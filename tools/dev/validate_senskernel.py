"""validate_senskernel.py - tailored sensitivity kernels (Swenson & Wahr 2002).

A basin average is  a = k' x  for some kernel k. The exact-indicator
kernel (shLowLevel.basinKernel) recovers the region perfectly IF the
data were perfect - but GRACE data are not, and the exact indicator has
energy at every degree, so it amplifies the high-degree noise without
limit. The standard remedy is to smooth, which trades noise for leakage
from outside the basin.

Swenson & Wahr (2002) make that trade explicitly rather than by picking
a filter radius. Minimise

    J(k) = (leakage)^2 + Alpha * (propagated noise)^2

    leakage^2 = || W (k - k_exact) ||^2   over the far field
    noise^2   = k' N k                    N = error covariance

which is a quadratic in k with the closed-form solution

    k = (M + Alpha N)^-1 M k_exact,   M = the far-field Gram matrix.

Alpha sweeps the trade-off: 0 gives the exact indicator (minimum
leakage, maximum noise), large Alpha gives a heavily smoothed kernel.
The GravIS/CCI kernels (Groh & Horwath 2021; Doehne et al. 2023) are
this construction with a particular weighting.

Checks: the two limits, that the trade-off curve is monotone in both
quantities (so Alpha really is a dial), and that an L-curve corner
exists and beats a Gaussian of matched noise.

Developed by Matthias Weigelt with the help of Claude (Opus 5),
2026-08-11 (v3.7.0).
"""
import numpy as np


def setup(lmax=20, seed=2):
    """Degree-only ("isotropic") toy: one coefficient per degree.

    Enough to show the trade-off honestly: the kernel is a vector of
    degree weights, leakage is its mismatch to the exact indicator
    outside the cap, and the noise grows with degree as GRACE's does.
    """
    n = np.arange(lmax + 1)
    # exact indicator of a spherical cap of radius psi0: Legendre coeffs
    psi0 = np.radians(10.0)
    x = np.cos(psi0)
    k_exact = np.zeros(lmax + 1)
    k_exact[0] = (1 - x) / 2
    Pnm1, Pn = np.ones_like(x), x
    P = [np.ones_like(x), x]
    for d in range(2, lmax + 2):
        Pn, Pnm1 = ((2 * d - 1) * x * Pn - (d - 1) * Pnm1) / d, Pn
        P.append(Pn)
    for d in range(1, lmax + 1):
        k_exact[d] = (P[d - 1] - P[d + 1]) / 2
    # GRACE-like error growth with degree
    rng = np.random.default_rng(seed)
    sig = 1e-11 * (1 + (n / 8.0) ** 3)
    N = np.diag(sig ** 2)
    # far-field weighting: degrees are all "outside" to some extent; use
    # the degree-amplitude of a unit far-field signal
    M = np.diag(1.0 / (1.0 + n))
    return n, k_exact, M, N, rng


def kernel(alpha, k_exact, M, N):
    """Minimise J subject to the UNIT-RESPONSE constraint k'kex = kex'kex.

    Without the constraint the cheapest way to cut noise is to shrink k
    towards zero, so the 'optimal' kernel measures nothing. One Lagrange
    multiplier restores it.
    """
    A = M + alpha * N
    k0 = np.linalg.solve(A, M @ k_exact)
    v = np.linalg.solve(A, k_exact)
    lam = (k_exact @ k_exact - k_exact @ k0) / (k_exact @ v)
    return k0 + lam * v


def metrics(k, k_exact, M, N):
    d = k - k_exact
    leak = np.sqrt(d @ M @ d)
    noise = np.sqrt(k @ N @ k)
    return leak, noise


def main():
    n, kex, M, N, rng = setup()
    print("exact indicator: leakage %.3e, noise %.3e"
          % metrics(kex, kex, M, N))

    print("\ntrade-off curve")
    alphas = np.logspace(14, 24, 11)
    L, S = [], []
    for a in alphas:
        k = kernel(a, kex, M, N)
        l, s = metrics(k, kex, M, N)
        L.append(l); S.append(s)
        print("  alpha %.1e: leakage %.3e  noise %.3e" % (a, l, s))
    L, S = np.array(L), np.array(S)

    # Alpha must be a genuine dial: leakage up, noise down, monotonically
    assert np.all(np.diff(L) > 0), "leakage must grow with alpha"
    assert np.all(np.diff(S) < 0), "noise must fall with alpha"
    print("\nmonotone in both: alpha is a real trade-off parameter")

    # the two limits
    k0 = kernel(0.0, kex, M, N)
    assert np.allclose(k0, kex, atol=1e-12), "alpha=0 must give the indicator"
    for a in (1e14, 1e20, 1e26):
        g = kernel(a, kex, M, N) @ kex / (kex @ kex)
        assert abs(g - 1.0) < 1e-8, (a, g)
    print("limits: alpha=0 -> exact indicator; gain stays 1 at every alpha")

    # an L-curve corner exists (maximum curvature in log-log)
    lx, ly = np.log(L), np.log(S)
    d1x, d1y = np.gradient(lx), np.gradient(ly)
    d2x, d2y = np.gradient(d1x), np.gradient(d1y)
    curv = np.abs(d1x * d2y - d1y * d2x) / (d1x ** 2 + d1y ** 2) ** 1.5
    ic = int(np.argmax(curv[1:-1])) + 1
    print("\nL-curve corner at alpha %.1e: leakage %.3e noise %.3e"
          % (alphas[ic], L[ic], S[ic]))
    assert 0 < ic < len(alphas) - 1, "the corner must be interior"

    # compare against a Gaussian kernel with the SAME noise: the tailored
    # kernel must leak less, which is the whole claim of the method
    target = S[ic]
    best = None
    for r in np.linspace(50, 3000, 200):
        b = np.log(2) / (1 - np.cos(r / 6371.0))
        w = np.exp(-n * (n + 1) / (2 * b))       # smooth Gaussian-like taper
        kg = kex * w
        kg = kg * (kex @ kex) / (kg @ kex)      # match the unit gain too
        lg, sg = metrics(kg, kex, M, N)
        if best is None or abs(sg - target) < abs(best[2] - target):
            best = (r, lg, sg)
    print("Gaussian at matched noise (r = %.0f km): leakage %.3e "
          "vs tailored %.3e" % (best[0], best[1], L[ic]))
    assert L[ic] < best[1], "the tailored kernel must beat a Gaussian"
    print("  tailored leaks %.1f%% less at the same noise"
          % (100 * (1 - L[ic] / best[1])))
    print("\n  The margin depends strongly on the far-field weighting:")
    for label, Mv in (("1/(n+1)^2", np.diag(1.0 / (n + 1.0) ** 2)),
                      ("1/(n+1)", M),
                      ("1", np.diag(np.ones_like(n, dtype=float))),
                      ("n+1", np.diag(n + 1.0))):
        Lw, Sw = [], []
        for a2 in alphas:
            kk = kernel(a2, kex, Mv, N)
            lw, sw = metrics(kk, kex, Mv, N)
            Lw.append(lw); Sw.append(sw)
        Lw, Sw = np.array(Lw), np.array(Sw)
        lxw, lyw = np.log(Lw), np.log(Sw)
        g1x, g1y = np.gradient(lxw), np.gradient(lyw)
        g2x, g2y = np.gradient(g1x), np.gradient(g1y)
        cw = np.abs(g1x * g2y - g1y * g2x) / (g1x ** 2 + g1y ** 2) ** 1.5
        iw = int(np.argmax(cw[1:-1])) + 1
        bb = None
        for r in np.linspace(50, 4000, 400):
            bq = np.log(2) / (1 - np.cos(r / 6371.0))
            w2 = np.exp(-n * (n + 1) / (2 * bq))
            kg2 = kex * w2
            kg2 = kg2 * (kex @ kex) / (kg2 @ kex)
            l2, s2 = metrics(kg2, kex, Mv, N)
            if bb is None or abs(s2 - Sw[iw]) < abs(bb[1] - Sw[iw]):
                bb = (l2, s2)
        print("    M = %-10s tailored %.3e vs gauss %.3e -> %.1f%% less"
              % (label, Lw[iw], bb[0], 100 * (1 - Lw[iw] / bb[0])))
    print("  The method is not magic: 2-16% here, and the ordering is "
          "not the intuitive one - the largest margin comes with the "
          "HIGH-degree weighting M = n+1, where a Gaussian's fixed "
          "shape is furthest from optimal. Report the number for YOUR "
          "weighting, measured, not a headline.")

    print("\nvalidate_senskernel: trade-off monotone, limits correct, "
          "beats a Gaussian at matched noise")


if __name__ == "__main__":
    main()
