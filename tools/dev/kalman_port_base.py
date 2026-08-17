"""Importable core of validate_kalman.py (shared by validate_kalman_qc.py).
Auto-derived; edit validate_kalman.py first. Claude (Fable 5), 2026-08-17."""
import numpy as np


# ---------------------------------------------------------------- VAR(p)
def empirical_cov(X, h):
    """Unbiased empirical (cross-)covariance Sigma(h) of centered series.

    X : (P, T) state realization, already centered/detrended.
    h : lag >= 0.  Sigma(h) = 1/(T-h) * sum_i x_i x_{i-h}^T
    """
    P, T = X.shape
    if h == 0:
        return X @ X.T / T
    return X[:, h:] @ X[:, :-h].T / (T - h)


def estimate_var(X, p=1, shrink=0.0):
    """VAR(p) by Yule-Walker.  Returns (Phi list, Q, Sigma0).

    Solves  [Sigma(1) ... Sigma(p)] = [Phi_1 ... Phi_p] * S
    with S the block-Toeplitz matrix S[i,j] = Sigma(i-j) (Sigma(-h)=Sigma(h)^T),
    then  Q = Sigma(0) - sum_k Phi_k Sigma(k)^T   (Luetkepohl 2005).
    p=1: Phi_1 = Sigma(1) Sigma(0)^-1, Q = Sigma(0) - Phi_1 Sigma(0) Phi_1^T
    == Kurtenbach (3.84)-(3.85) with Sigma_bar = Sigma(0), Sigma_Delta = Sigma(1).
    shrink : diagonal loading factor on Sigma(0) blocks (robustness).
    """
    P = X.shape[0]
    Sig = [empirical_cov(X, h) for h in range(p + 1)]
    if shrink > 0:
        Sig[0] = Sig[0] + shrink * np.trace(Sig[0]) / P * np.eye(P)
    # block-Toeplitz S (p*P x p*P)
    S = np.empty((p * P, p * P))
    for i in range(p):
        for j in range(p):
            h = i - j
            S[i*P:(i+1)*P, j*P:(j+1)*P] = Sig[h] if h >= 0 else Sig[-h].T
    G = np.hstack([Sig[h] for h in range(1, p + 1)])          # (P, p*P)
    Phi_all = np.linalg.solve(S.T, G.T).T                      # (P, p*P)
    Phi = [Phi_all[:, k*P:(k+1)*P] for k in range(p)]
    Q = Sig[0] - sum(Phi[k] @ Sig[k+1].T for k in range(p))
    Q = 0.5 * (Q + Q.T)
    return Phi, Q, Sig[0]


def companion(Phi, Q):
    """VAR(p) -> VAR(1) companion form (Kvas Sec. 2.3.1)."""
    p = len(Phi)
    P = Phi[0].shape[0]
    B = np.zeros((p * P, p * P))
    B[:P, :] = np.hstack(Phi)
    if p > 1:
        B[P:, :-P] = np.eye((p - 1) * P)
    Qt = np.zeros((p * P, p * P))
    Qt[:P, :P] = Q
    return B, Qt


def stationary_cov(B, Qt, n_iter=2000, tol=1e-13):
    """Solve discrete Lyapunov S = B S B' + Qt by fixed-point iteration."""
    S = Qt.copy()
    for _ in range(n_iter):
        Sn = B @ S @ B.T + Qt
        if np.max(np.abs(Sn - S)) < tol * max(1.0, np.max(np.abs(S))):
            return 0.5 * (Sn + Sn.T)
        S = Sn
    return 0.5 * (S + S.T)


# ------------------------------------------------------- Kalman filter/RTS
def kalman_filter(B, Q, Sigma0, obs):
    """Forward Kalman filter with stationary initialization.

    obs : list over epochs; each entry is
          None                        -> prediction only (data gap)
          ('sol', l, R)               -> l = x + v, v ~ N(0, R)
          ('neq', N, b)               -> normal equation N x = b (information update)
    Returns dict with xp/Pp (predicted), xf/Pf (filtered), contrib
    (diagonal share of data in the estimate: diag(K) resp. diag(Pf N)).
    """
    P = B.shape[0]
    T = len(obs)
    xp = np.zeros((P, T)); xf = np.zeros((P, T))
    Pp = np.zeros((T, P, P)); Pf = np.zeros((T, P, P))
    contrib = np.full((P, T), np.nan)
    x_prev = np.zeros(P)
    P_prev = Sigma0.copy()       # so that Pp[0] = Q + B Sigma0 B' = Sigma0
    for t in range(T):
        # predict (Kvas 2.90-2.92: with P_{-1}=Sigma0 the prediction is Sigma0)
        x_m = B @ x_prev
        P_m = B @ P_prev @ B.T + Q
        P_m = 0.5 * (P_m + P_m.T)
        xp[:, t] = x_m; Pp[t] = P_m
        o = obs[t]
        if o is None:
            x_pl, P_pl = x_m, P_m
        elif o[0] == 'sol':
            _, l, R = o
            K = np.linalg.solve((P_m + R).T, P_m.T).T          # P_m (P_m+R)^-1
            x_pl = x_m + K @ (l - x_m)
            P_pl = (np.eye(P) - K) @ P_m
            contrib[:, t] = np.diag(K)
        elif o[0] == 'neq':
            _, N, b = o
            Jm = np.linalg.inv(P_m)                            # information form
            P_pl = np.linalg.inv(Jm + N)
            x_pl = P_pl @ (Jm @ x_m + b)
            contrib[:, t] = np.diag(P_pl @ N)
        else:
            raise ValueError(o[0])
        P_pl = 0.5 * (P_pl + P_pl.T)
        xf[:, t] = x_pl; Pf[t] = P_pl
        x_prev, P_prev = x_pl, P_pl
    return dict(xp=xp, Pp=Pp, xf=xf, Pf=Pf, contrib=contrib)


def rts_smoother(B, filt):
    """RTS backward pass: G_t = Pf_t B' Pp_{t+1}^-1."""
    xf, Pf, xp, Pp = filt['xf'], filt['Pf'], filt['xp'], filt['Pp']
    P, T = xf.shape
    xs = xf.copy(); Ps = Pf.copy()
    for t in range(T - 2, -1, -1):
        G = np.linalg.solve(Pp[t+1].T, (Pf[t] @ B.T).T).T
        xs[:, t] = xf[:, t] + G @ (xs[:, t+1] - xp[:, t+1])
        Ps[t] = Pf[t] + G @ (Ps[t+1] - Pp[t+1]) @ G.T
        Ps[t] = 0.5 * (Ps[t] + Ps[t].T)
    return dict(xs=xs, Ps=Ps)


def batch_lsa(B, Q, Sigma0, obs):
    """Joint least squares over all epochs (Kvas Sec. 2.3): the reference.

    Pseudo-observations: x_0 ~ N(0, Sigma0); x_t - B x_{t-1} ~ N(0, Q);
    plus the per-epoch data.  Normal matrix is block tridiagonal.
    """
    P = B.shape[0]; T = len(obs)
    N = np.zeros((T * P, T * P)); rhs = np.zeros(T * P)
    iS0 = np.linalg.inv(Sigma0); iQ = np.linalg.inv(Q)
    N[:P, :P] += iS0
    for t in range(1, T):
        i0, i1 = (t-1)*P, t*P
        N[i1:i1+P, i1:i1+P] += iQ
        N[i0:i0+P, i0:i0+P] += B.T @ iQ @ B
        N[i0:i0+P, i1:i1+P] += -B.T @ iQ
        N[i1:i1+P, i0:i0+P] += -iQ @ B
    for t, o in enumerate(obs):
        if o is None:
            continue
        i0 = t * P
        if o[0] == 'sol':
            _, l, R = o
            iR = np.linalg.inv(R)
            N[i0:i0+P, i0:i0+P] += iR
            rhs[i0:i0+P] += iR @ l
        else:
            _, Nt, bt = o
            N[i0:i0+P, i0:i0+P] += Nt
            rhs[i0:i0+P] += bt
    x = np.linalg.solve(N, rhs)
    return x.reshape(T, P).T


# ------------------------------------------------------------- taper
def taper_cov(C, psi, psi0):
    """Kvas eq. (2.120): scale correlations by exp(-psi/psi0).

    C   : (M, M) spatial covariance
    psi : (M, M) spherical distances [rad or km, consistent with psi0]
    """
    if np.isinf(psi0):
        return C.copy()
    return C * np.exp(-psi / psi0)


