"""Python validation port for shAnalysis Kalman/VAR module.

Implements and cross-validates, before any MATLAB code exists:
  - estimate_var : VAR(p) via Yule-Walker (Kvas 2019, Sec. 2.4);
                   p=1 reduces to Kurtenbach 2012 eqs. (3.84)-(3.88)
  - kalman_filter: forward filter, 'solution' and 'neq' (information) modes,
                   stationary init x0-=0, P0-=Sigma(0) (Kvas eqs. 2.90-2.92)
  - rts_smoother : Rauch-Tung-Striebel backward pass
  - batch_lsa    : block-tridiagonal joint least squares (Kvas Sec. 2.3),
                   the reference the smoother must reproduce exactly
  - taper_cov    : distance-dependent correlation scaling (Kvas eq. 2.120)

Prepared by Claude (Fable 5), 2026-08-17.
"""
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


# ================================================================= tests
def _sim_var(Phi, Q, T, rng, burn=500):
    p = len(Phi); P = Phi[0].shape[0]
    L = np.linalg.cholesky(Q + 1e-15 * np.eye(P))
    X = np.zeros((P, T + burn))
    for t in range(p, T + burn):
        X[:, t] = sum(Phi[k] @ X[:, t-1-k] for k in range(p)) + L @ rng.standard_normal(P)
    return X[:, burn:]


def run_tests():
    rng = np.random.default_rng(42)
    P, T = 6, 4000
    A = rng.standard_normal((P, P))
    Phi_true = [0.9 * A / np.max(np.abs(np.linalg.eigvals(A)))]
    Qt_true = np.eye(P) * 0.2
    X = _sim_var(Phi_true, Qt_true, T, rng)

    # T1: VAR(1) Yule-Walker == Kurtenbach closed form (3.84)-(3.85)
    Phi, Q, S0 = estimate_var(X, 1)
    S1 = empirical_cov(X, 1)
    B_kurt = np.linalg.solve(S0.T, S1.T).T
    Q_kurt = S0 - B_kurt @ S0 @ B_kurt.T
    assert np.allclose(Phi[0], B_kurt, atol=1e-12), "T1 Phi"
    assert np.allclose(Q, 0.5*(Q_kurt+Q_kurt.T), atol=1e-12), "T1 Q"
    print("T1  VAR(1) YW == Kurtenbach closed form            PASS")

    # T2: Q symmetric PSD, companion spectral radius < 1
    B, Qc = companion(Phi, Q)
    ev = np.linalg.eigvalsh(Q)
    assert ev.min() > -1e-10 * ev.max(), "T2 PSD"
    assert np.max(np.abs(np.linalg.eigvals(B))) < 1, "T2 stability"
    print("T2  Q PSD, companion stable                        PASS")

    # T3: KF + RTS == batch LSA (Kvas Sec. 2.3), with a data gap
    Sig_st = stationary_cov(B, Qc)
    Tobs = 12
    obs = []
    for t in range(Tobs):
        if t == 5:
            obs.append(None)                                   # gap
        else:
            R = np.diag(0.3 + rng.random(P))
            l = rng.standard_normal(P)
            obs.append(('sol', l, R))
    filt = kalman_filter(B, Qc, Sig_st, obs)
    smo = rts_smoother(B, filt)
    xb = batch_lsa(B, Qc, Sig_st, obs)
    err = np.max(np.abs(smo['xs'] - xb)) / np.max(np.abs(xb))
    assert err < 1e-9, f"T3 rel err {err:.2e}"
    print(f"T3  KF+RTS == batch LSA (rel {err:.1e}, incl. gap) PASS")

    # T4: NEQ (information) mode == solution mode when N=R^-1, b=N l
    obs_neq = []
    for o in obs:
        if o is None:
            obs_neq.append(None)
        else:
            _, l, R = o
            Nt = np.linalg.inv(R)
            obs_neq.append(('neq', Nt, Nt @ l))
    f2 = kalman_filter(B, Qc, Sig_st, obs_neq)
    assert np.allclose(f2['xf'], filt['xf'], atol=1e-9), "T4 states"
    assert np.allclose(f2['Pf'], filt['Pf'], atol=1e-9), "T4 covs"
    print("T4  NEQ mode == solution mode                      PASS")

    # T5: VAR(2) coefficient recovery from long simulation
    P2 = 4
    Phi2_true = [np.diag(rng.uniform(0.3, 0.5, P2)),
                 np.diag(rng.uniform(0.1, 0.3, P2))]
    Q2_true = np.eye(P2) * 0.1
    X2 = _sim_var(Phi2_true, Q2_true, 200000, rng)
    Phi2, Q2, _ = estimate_var(X2, 2)
    for k in range(2):
        assert np.max(np.abs(Phi2[k] - Phi2_true[k])) < 0.03, f"T5 Phi_{k+1}"
    assert np.max(np.abs(Q2 - Q2_true)) < 0.01, "T5 Q"
    print("T5  VAR(2) recovery via companion YW               PASS")

    # T6: gap -> filter covariance relaxes toward Sigma(0)
    obs_gap = [obs[0]] + [None] * 60
    fg = kalman_filter(B, Qc, Sig_st, obs_gap)
    d_end = np.max(np.abs(fg['Pf'][-1] - Sig_st)) / np.max(np.abs(Sig_st))
    d_start = np.max(np.abs(fg['Pf'][0] - Sig_st)) / np.max(np.abs(Sig_st))
    assert d_end < 1e-3 and d_start > 0.05, "T6"
    print(f"T6  long gap relaxes to Sigma(0) ({d_start:.2f}->{d_end:.1e}) PASS")

    # T7: B=0 update == per-epoch Wiener filter (tvANS connection)
    B0 = np.zeros_like(Sig_st); Q0 = Sig_st
    f0 = kalman_filter(B0, Q0, Sig_st, obs[:1])
    _, l, R = obs[0]
    W = np.linalg.solve((Sig_st + R).T, Sig_st.T).T
    assert np.allclose(f0['xf'][:, 0], W @ l, atol=1e-10), "T7"
    print("T7  B=0 update == Wiener (tvANS consistency)       PASS")

    # T8: taper: psi0=inf identity; result symmetric PSD for exp taper
    M = 20
    pts = rng.standard_normal((3, M)); pts /= np.linalg.norm(pts, axis=0)
    psi = np.arccos(np.clip(pts.T @ pts, -1, 1))
    Craw = empirical_cov(rng.standard_normal((M, 500)), 0)
    assert np.allclose(taper_cov(Craw, psi, np.inf), Craw), "T8 identity"
    Ct = taper_cov(Craw, psi, 0.5)
    assert np.allclose(Ct, Ct.T), "T8 sym"
    # exp(-psi/psi0) is a PSD kernel on the sphere -> Schur product PSD
    assert np.linalg.eigvalsh(Ct).min() > -1e-8 * np.abs(Ct).max(), "T8 PSD"
    print("T8  taper: identity at inf, Schur-product PSD      PASS")

    print("\nAll 8 tests PASS")


if __name__ == "__main__":
    run_tests()
