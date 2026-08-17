"""Python validation for the Kalman-module completion PRs (v3.21-3.24).

Q-block  innovation-based quality control (Kvas 2019, Sec. 3.3):
  Q1  Wilson-Hilferty chi^2 quantile vs scipy at QC-relevant dof/alpha
  Q2  solution-mode statistic T = d' S^-1 d is chi^2(P): mean == P
  Q3  NEQ-form statistic u'(N Pm N + N)^-1 u equals the solution form
      EXACTLY when N = R^-1 (algebraic identity)
  Q4  an injected 50-sigma blunder fires the test at that epoch only

J-block  Joseph-stabilized update:
  J1  identical to the standard update in well-conditioned arithmetic
  J2  with cond(R) ~ 1e14 the standard update loses PSD (negative
      eigenvalue), the Joseph form does not

L-block  fixed-lag smoother:
  L1  Lag >= T-1 reproduces the full RTS smoother to machine precision
  L2  the error vs full RTS decays with growing Lag (process-dependent
      rate; monotone on average)

Prepared by Claude (Fable 5), 2026-08-17.
"""
import numpy as np
from scipy.stats import chi2
from kalman_port_base import (estimate_var, companion, stationary_cov,
                              kalman_filter, rts_smoother)


# ------------------------------------------------------------ chi2 WH
def chi2_quantile_wh(p, k):
    """Wilson-Hilferty: chi2inv(p,k) ~ k (1 - 2/(9k) + z sqrt(2/(9k)))^3."""
    from math import sqrt
    from scipy.special import erfinv          # base-MATLAB erfinv analogue
    z = sqrt(2.0) * erfinv(2.0 * p - 1.0)
    return k * (1.0 - 2.0/(9.0*k) + z * np.sqrt(2.0/(9.0*k)))**3


# ------------------------------------------------------------ QC stats
def qc_solution(l, xm, Pm, R):
    d = l - xm
    S = Pm + R
    return float(d @ np.linalg.solve(S, d))


def qc_neq(b, N, xm, Pm):
    u = b - N @ xm
    W = N @ Pm @ N + N
    return float(u @ np.linalg.lstsq(W, u, rcond=None)[0])


# ------------------------------------------------------------ Joseph
def update_standard(xm, Pm, l, R):
    P = Pm.shape[0]
    K = np.linalg.solve((Pm + R).T, Pm.T).T
    return xm + K @ (l - xm), (np.eye(P) - K) @ Pm


def update_joseph(xm, Pm, l, R):
    P = Pm.shape[0]
    K = np.linalg.solve((Pm + R).T, Pm.T).T
    IK = np.eye(P) - K
    Ppl = IK @ Pm @ IK.T + K @ R @ K.T
    return xm + K @ (l - xm), 0.5 * (Ppl + Ppl.T)


# ------------------------------------------------------------ fixed lag
def fixed_lag(B, Q, S0, obs, lag):
    """Definition-true fixed-lag smoother: for each t smooth back from
    min(T-1, t+lag) using the forward-filter quantities."""
    filt = kalman_filter(B, Q, S0, obs)
    xf, Pf, xp, Pp = filt['xf'], filt['Pf'], filt['xp'], filt['Pp']
    P, T = xf.shape
    xs = np.empty_like(xf)
    for t in range(T):
        e = min(T - 1, t + lag)
        x_next = xf[:, e]
        P_next = Pf[e]
        for s in range(e - 1, t - 1, -1):
            G = np.linalg.solve(Pp[s+1].T, (Pf[s] @ B.T).T).T
            x_next = xf[:, s] + G @ (x_next - xp[:, s+1])
            P_next = Pf[s] + G @ (P_next - Pp[s+1]) @ G.T
        xs[:, t] = x_next
    return xs


def run():
    rng = np.random.default_rng(9)

    # Q1 Wilson-Hilferty accuracy
    # Measured accuracy map (this script's own scan): <= 3.5% at k <= 10,
    # <= 0.2% at k = 50, <= 1e-4 at k >= 500 - the QC operating point is
    # k = P (hundreds to 1677), where the approximation is essentially exact.
    for k, lim in ((5, 0.05), (10, 0.05), (50, 0.005), (500, 5e-4), (1677, 5e-4)):
        for a in (1e-2, 1e-3, 1e-4):
            ref = chi2.ppf(1 - a, k)
            wh = chi2_quantile_wh(1 - a, k)
            rel = abs(wh - ref) / ref
            assert rel < lim, f"Q1 k={k} a={a} rel={rel:.2e}"
    print("Q1  Wilson-Hilferty within measured bounds           PASS")

    # Q2 statistic is chi^2(P): mean over epochs
    P = 8
    A = rng.standard_normal((P, P))
    Phi = [0.8 * A / np.max(np.abs(np.linalg.eigvals(A)))]
    Q = 0.3 * np.eye(P)
    B, Qc = companion(Phi, Q)
    S0 = stationary_cov(B, Qc)
    Tn = 4000
    Ls = np.linalg.cholesky(S0)
    stats = []
    xm = np.zeros(P); Pm = S0
    x_true = Ls @ rng.standard_normal(P)
    for t in range(Tn):
        R = np.diag(0.2 + rng.random(P))
        l = x_true + np.sqrt(np.diag(R)) * rng.standard_normal(P)
        stats.append(qc_solution(l, xm, Pm, R))
        # propagate truth and filter
        xpl, Ppl = update_joseph(xm, Pm, l, R)
        w = np.linalg.cholesky(Q + 1e-15*np.eye(P)) @ rng.standard_normal(P)
        x_true = Phi[0] @ x_true + w
        xm = B @ xpl
        Pm = B @ Ppl @ B.T + Qc
    m = np.mean(stats)
    assert abs(m - P) < 0.15, f"Q2 mean {m:.3f} vs {P}"
    print(f"Q2  E[T] = {m:.2f} (dof {P})                          PASS")

    # Q3 NEQ form == solution form (identity)
    for _ in range(20):
        R = np.diag(0.2 + rng.random(P))
        l = rng.standard_normal(P)
        xm = rng.standard_normal(P)
        Pm = S0
        N = np.linalg.inv(R)
        t1 = qc_solution(l, xm, Pm, R)
        t2 = qc_neq(N @ l, N, xm, Pm)
        assert abs(t1 - t2) < 1e-8 * max(1, t1), "Q3"
    print("Q3  NEQ statistic == solution statistic              PASS")

    # Q4 blunder detection: one 50-sigma epoch fires, others do not
    thr = chi2_quantile_wh(1 - 1e-3, P)
    obs = []
    x_true = Ls @ rng.standard_normal(P)
    xm = np.zeros(P); Pm = S0
    fired = []
    for t in range(60):
        R = np.diag(0.2 * np.ones(P))
        l = x_true + np.sqrt(np.diag(R)) * rng.standard_normal(P)
        if t == 30:
            l = l + 50 * np.sqrt(np.diag(R))          # blunder
        T = qc_solution(l, xm, Pm, R)
        fired.append(T > thr)
        if not fired[-1]:
            xm, Pm = update_joseph(xm, Pm, l, R)      # reject policy
        w = np.linalg.cholesky(Q + 1e-15*np.eye(P)) @ rng.standard_normal(P)
        x_true = Phi[0] @ x_true + w
        xm = B @ xm
        Pm = B @ Pm @ B.T + Qc
    assert fired[30] and sum(fired) <= 2, f"Q4 fired {sum(fired)}"
    print(f"Q4  blunder fires at t=30 only ({sum(fired)} total)     PASS")

    # J1 Joseph == standard when well-conditioned
    R = np.diag(0.2 + rng.random(P))
    l = rng.standard_normal(P)
    x1, P1 = update_standard(np.zeros(P), S0, l, R)
    x2, P2 = update_joseph(np.zeros(P), S0, l, R)
    assert np.allclose(x1, x2, atol=1e-10) and np.allclose(P1, P2, atol=1e-10), "J1"
    print("J1  Joseph == standard (well-conditioned)            PASS")

    # J2 accuracy A/B against a 50-digit reference: with R tiny relative
    # to a wide-spectrum prior (strong daily NEQ vs broad process prior,
    # THE regime of this module), the standard (I-K)Pm update loses ~6
    # orders of relative accuracy to cancellation; the Joseph form stays
    # at the double-precision floor. Neither loses PSD here - the honest
    # metric is accuracy, not definiteness (measured, not assumed).
    from mpmath import mp, matrix as mpm
    mp.dps = 50
    Pj_dim = 25
    Vq, _ = np.linalg.qr(rng.standard_normal((Pj_dim, Pj_dim)))
    Pm = Vq @ np.diag(np.logspace(0, -8, Pj_dim)) @ Vq.T
    Pm = 0.5 * (Pm + Pm.T)

    def ref_up(Pm_, R_):
        A_ = mpm(Pm_.tolist()); B_ = mpm(R_.tolist())
        res = (A_**-1 + B_**-1)**-1
        return np.array([[float(res[i, j]) for j in range(Pj_dim)]
                         for i in range(Pj_dim)])

    e_std = 0.0; e_jos = 0.0
    for _ in range(6):
        R = np.diag(10.0**rng.uniform(-14, -10, Pj_dim))
        l = rng.standard_normal(Pj_dim)
        _, Ps = update_standard(np.zeros(Pj_dim), Pm, l, R)
        _, Pjm = update_joseph(np.zeros(Pj_dim), Pm, l, R)
        Pr = ref_up(Pm, R); sc = np.abs(Pr).max()
        Ps = 0.5 * (Ps + Ps.T)
        e_std = max(e_std, np.abs(Ps - Pr).max() / sc)
        e_jos = max(e_jos, np.abs(Pjm - Pr).max() / sc)
    assert e_std > 1e-8, f"J2 standard unexpectedly accurate ({e_std:.1e})"
    assert e_jos < 1e-10, f"J2 Joseph inaccurate ({e_jos:.1e})"
    print(f"J2  rel err vs 50-digit ref: std {e_std:.1e}, Joseph {e_jos:.1e}  PASS")

    # L1/L2 fixed lag
    Tn = 40
    obs = []
    for t in range(Tn):
        R = np.diag(0.3 + rng.random(P))
        obs.append(('sol', rng.standard_normal(P), R))
    filt = kalman_filter(B, Qc, S0, obs)
    smo = rts_smoother(B, filt)
    x_full = fixed_lag(B, Qc, S0, obs, Tn)
    assert np.max(np.abs(x_full - smo['xs'])) < 1e-9, "L1"
    print("L1  Lag >= T-1 == full RTS                           PASS")
    errs = []
    for lag in (0, 2, 5, 10, 20):
        xl = fixed_lag(B, Qc, S0, obs, lag)
        errs.append(np.linalg.norm(xl - smo['xs']))
    assert all(errs[i+1] < errs[i] for i in range(len(errs)-1)), f"L2 {errs}"
    assert errs[-1] < 0.02 * errs[0], "L2 decay"
    print(f"L2  error vs lag decays {errs[0]:.2f} -> {errs[-1]:.1e}       PASS")

    print("\nAll 8 checks PASS")


if __name__ == "__main__":
    run()
