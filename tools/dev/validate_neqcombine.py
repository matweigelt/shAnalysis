"""Python validation for shLowLevel.neqCombine (NEQ-level VCE).

Combination of K normal equations on the normal-equation level (the
COST-G approach; cf. Kvas thesis Sec. 2.1-2.2 for the NEQ formalism)
with one variance component per contribution, Foerstner/Koch iteration:

    N      = sum_i (1/s_i^2) N_i,   b = sum_i (1/s_i^2) b_i
    x      = N^-1 b
    Om_i   = ltpl_i - 2 x' b_i + x' N_i x        (= ||l_i - A_i x||^2)
    r_i    = n_i - (1/s_i^2) tr(N^-1 N_i)        (partial redundancy)
    s_i^2 <- Om_i / r_i

Checks:
  V1  fixed weights: combined x equals the directly stacked GLS solve
  V2  VCE recovers known variance factors [1, 4] within scatter
  V3  redundancy invariant sum r_i = sum n_i - P at every iteration
  V4  combined solution beats every single contribution against truth
  V5  identical noise -> variance factors agree

Prepared by Claude (Fable 5), 2026-08-17.
"""
import numpy as np


def neq_combine(neqs, weights=None, max_iter=20, tol=1e-6):
    """neqs: list of dicts N, b, ltpl, nobs. Returns x, info."""
    K = len(neqs)
    P = neqs[0]['N'].shape[0]
    if weights is not None:
        s2 = 1.0 / np.asarray(weights, float)
        fixed = True
    else:
        s2 = np.ones(K)
        fixed = False
    it = 0
    conv = fixed
    red_hist = []
    for it in range(1, max_iter + 1):
        N = sum(neqs[i]['N'] / s2[i] for i in range(K))
        b = sum(neqs[i]['b'] / s2[i] for i in range(K))
        x = np.linalg.solve(N, b)
        if fixed:
            break
        Ni_inv = np.linalg.inv(N)
        s2_new = np.empty(K)
        r = np.empty(K)
        for i in range(K):
            Om = neqs[i]['ltpl'] - 2 * x @ neqs[i]['b'] + x @ neqs[i]['N'] @ x
            r[i] = neqs[i]['nobs'] - np.trace(Ni_inv @ neqs[i]['N']) / s2[i]
            s2_new[i] = Om / r[i]
        red_hist.append(r.copy())
        if np.max(np.abs(s2_new - s2) / s2) < tol:
            s2 = s2_new
            conv = True
            break
        s2 = s2_new
    N = sum(neqs[i]['N'] / s2[i] for i in range(K))
    b = sum(neqs[i]['b'] / s2[i] for i in range(K))
    x = np.linalg.solve(N, b)
    return x, dict(sigma2=s2, nIter=it, converged=conv, N=N,
                   red_hist=red_hist)


def make_group(A, x_true, sigma, rng):
    n = A.shape[0]
    l = A @ x_true + sigma * rng.standard_normal(n)
    return dict(N=A.T @ A, b=A.T @ l, ltpl=float(l @ l), nobs=n,
                A=A, l=l)


def run():
    rng = np.random.default_rng(4)
    P = 12
    x_true = rng.standard_normal(P)
    sig_true = [1.0, 2.0]                       # variance factors 1 and 4
    groups = []
    for k, s in enumerate(sig_true):
        A = rng.standard_normal((4000, P))
        groups.append(make_group(A, x_true, s, rng))

    # V1 fixed weights == stacked GLS
    w = [1.0, 0.25]
    x_fix, _ = neq_combine(groups, weights=w)
    Astk = np.vstack([np.sqrt(w[i]) * groups[i]['A'] for i in range(2)])
    lstk = np.concatenate([np.sqrt(w[i]) * groups[i]['l'] for i in range(2)])
    x_ref = np.linalg.lstsq(Astk, lstk, rcond=None)[0]
    assert np.max(np.abs(x_fix - x_ref)) < 1e-9, "V1"
    print("V1  fixed weights == stacked GLS                 PASS")

    # V2 VCE recovers [1, 4]
    x_vce, info = neq_combine(groups)
    s2 = info['sigma2']
    assert info['converged'], "V2 conv"
    assert abs(s2[0] - 1.0) < 0.08 and abs(s2[1] - 4.0) < 0.3, f"V2 {s2}"
    print(f"V2  VCE factors {s2[0]:.3f}, {s2[1]:.3f} (true 1, 4)      PASS")

    # V3 redundancy invariant at every iteration
    n_tot = sum(g['nobs'] for g in groups)
    for r in info['red_hist']:
        assert abs(r.sum() - (n_tot - P)) < 1e-6 * n_tot, "V3"
    print("V3  sum r_i = sum n_i - P at every iteration     PASS")

    # V4 combined beats each single contribution
    e_comb = np.linalg.norm(x_vce - x_true)
    for g in groups:
        xi = np.linalg.solve(g['N'], g['b'])
        assert e_comb < np.linalg.norm(xi - x_true), "V4"
    print("V4  combined beats every single contribution     PASS")

    # V5 identical noise -> equal factors
    g2 = [make_group(rng.standard_normal((3000, P)), x_true, 1.5, rng)
          for _ in range(3)]
    _, info2 = neq_combine(g2)
    s2 = info2['sigma2']
    assert np.max(np.abs(s2 / s2.mean() - 1)) < 0.1, f"V5 {s2}"
    print("V5  identical noise -> equal variance factors    PASS")

    print("\nAll 5 checks PASS")


if __name__ == "__main__":
    run()
