"""validate_leakage.py - the leakage numerics, in Python, BEFORE MATLAB.

Two corrections for the signal a GRACE filter removes:

  forward modelling  iteratively find the mass field that, pushed through
                     the SAME filter as the data, reproduces the observed
                     filtered field. Fixed-point iteration
                         m_(k+1) = m_k + g * (obs - F(m_k))
                     with gain g. F is the full chain synthesis(filter(
                     analysis(.))), so it carries whatever filter the
                     data saw. Optionally the update is confined to a
                     region mask (mass is known to live only there).

  gridded scaling    per-pixel k(lat,lon) from a model series pushed
                     through the same chain:
                         k = sum_t(true*filt) / sum_t(filt^2)
                     the least-squares gain per pixel (Landerer/Swenson
                     style). Applied as corrected = k .* filtered.

Everything is checked against known truth on a synthetic cap, plus the
identities the MATLAB port must reproduce. Numbers printed here become
the pinned expectations in tests/testCorrectness.m.

Run: python3 tools/dev/validate_leakage.py

Developed by Matthias Weigelt with the help of Claude (Opus 5),
2026-08-11 (v3.2.0).
"""
import numpy as np

# --------------------------------------------------------------- SH core
def legendre_alf(nmax, theta):
    """Fully normalized associated Legendre functions P[n, m, i].

    Standard forward column recursion; theta in RADIANS (colatitude).
    """
    t = np.atleast_1d(theta).astype(float)
    c, s = np.cos(t), np.sin(t)
    P = np.zeros((nmax + 1, nmax + 1, t.size))
    P[0, 0] = 1.0
    if nmax >= 1:
        P[1, 0] = np.sqrt(3.0) * c
        P[1, 1] = np.sqrt(3.0) * s
    for n in range(2, nmax + 1):
        # sectorial
        P[n, n] = np.sqrt((2.0 * n + 1.0) / (2.0 * n)) * s * P[n - 1, n - 1]
        for m in range(0, n):
            a = np.sqrt((2.0 * n - 1.0) * (2.0 * n + 1.0)
                        / ((n - m) * (n + m)))
            b = np.sqrt((2.0 * n + 1.0) * (n + m - 1.0) * (n - m - 1.0)
                        / ((n - m) * (n + m) * (2.0 * n - 3.0)))
            P[n, m] = a * c * P[n - 1, m] - (b * P[n - 2, m] if n - 2 >= m
                                             else 0.0)
    return P


def gauss_grid(nmax):
    """Gauss-Legendre latitudes and longitudes, exact to degree nmax."""
    nlat = nmax + 1
    x, w = np.polynomial.legendre.leggauss(nlat)
    theta = np.arccos(x)                       # colatitude [rad]
    nlon = 2 * nmax + 2
    lam = 2.0 * np.pi * np.arange(nlon) / nlon
    return theta, w, lam


def synthesis(C, S, theta, lam):
    """Grid (nlat x nlon) from fully normalized coefficients."""
    nmax = C.shape[0] - 1
    P = legendre_alf(nmax, theta)
    m = np.arange(nmax + 1)
    cosm = np.cos(np.outer(m, lam))            # (nmax+1, nlon)
    sinm = np.sin(np.outer(m, lam))
    A = np.einsum("nmi,nm->mi", P, C)          # (nmax+1, nlat)
    B = np.einsum("nmi,nm->mi", P, S)
    return A.T @ cosm + B.T @ sinm             # (nlat, nlon)


def analysis(grid, theta, w, lam, nmax):
    """Coefficients from a Gauss grid - exact inverse of synthesis."""
    nlon = lam.size
    m = np.arange(nmax + 1)
    cosm = np.cos(np.outer(m, lam))
    sinm = np.sin(np.outer(m, lam))
    # longitude quadrature
    a = grid @ cosm.T * (2.0 / nlon)           # (nlat, nmax+1)
    b = grid @ sinm.T * (2.0 / nlon)
    # C_nm = (1/4pi) * int f Ybar dOmega, discretised with the Gauss
    # weights in x = cos(theta) and dlambda = 2pi/nlon, giving the SAME
    # factor for every order: no special case at m = 0. (Checked against
    # int Pbar_10^2 sin = 2 and int Pbar_11^2 sin = 4.)
    P = legendre_alf(nmax, theta)
    C = np.einsum("nmi,i,im->nm", P, w, a) * 0.25
    S = np.einsum("nmi,i,im->nm", P, w, b) * 0.25
    tri = np.tril(np.ones((nmax + 1, nmax + 1), dtype=bool))
    return C * tri, S * tri


def gaussian_weights(nmax, radius_km, R_km=6378.1363):
    """Jekeli recursion; W[0] = 1."""
    b = np.log(2.0) / (1.0 - np.cos(radius_km / R_km))
    W = np.zeros(nmax + 1)
    W[0] = 1.0
    if nmax >= 1:
        W[1] = (1.0 + np.exp(-2.0 * b)) / (1.0 - np.exp(-2.0 * b)) - 1.0 / b
    for n in range(1, nmax):
        W[n + 1] = -(2.0 * n + 1.0) / b * W[n] + W[n - 1]
        if W[n + 1] < 0:                       # recursion goes unstable
            W[n + 1:] = 0.0
            break
    return W


# --------------------------------------------------- the chain and truth
NMAX = 40
THETA, WQ, LAM = gauss_grid(NMAX)
LATD = 90.0 - np.degrees(THETA)
LOND = np.degrees(LAM)
GW = gaussian_weights(NMAX, 500.0)


def smooth(field):
    """The processing chain the data saw: analysis -> Gauss -> synthesis."""
    C, S = analysis(field, THETA, WQ, LAM, NMAX)
    C = C * GW[:, None]
    S = S * GW[:, None]
    return synthesis(C, S, THETA, LAM)


def cap(lat0, lon0, radius_deg, amp):
    """A uniform disc of mass - the classic known-truth leakage target."""
    la = np.radians(LATD)[:, None]
    lo = np.radians(LOND)[None, :]
    psi = np.arccos(np.clip(
        np.sin(np.radians(lat0)) * np.sin(la)
        + np.cos(np.radians(lat0)) * np.cos(la)
        * np.cos(lo - np.radians(lon0)), -1, 1))
    return np.where(psi <= np.radians(radius_deg), amp, 0.0)


# ---------------------------------------------------------- sanity first
def check_transforms():
    rng = np.random.default_rng(7)
    C = np.tril(rng.normal(size=(NMAX + 1, NMAX + 1))) * 1e-3
    S = np.tril(rng.normal(size=(NMAX + 1, NMAX + 1))) * 1e-3
    S[:, 0] = 0.0
    g = synthesis(C, S, THETA, LAM)
    C2, S2 = analysis(g, THETA, WQ, LAM, NMAX)
    eC = np.max(np.abs(C2 - C))
    eS = np.max(np.abs(S2 - S))
    print("analysis(synthesis(x)) == x :  maxerr C %.2e  S %.2e" % (eC, eS))
    assert eC < 1e-12 and eS < 1e-12, "round trip broken"
    print("gaussian W[0], W[1], W[nmax] : %.6f %.6f %.6f"
          % (GW[0], GW[1], GW[NMAX]))
    assert abs(GW[0] - 1.0) < 1e-14 and np.all(np.diff(GW) <= 1e-15)


# -------------------------------------------------- forward modelling
def forward_model(obs, mask=None, gain=1.0, max_iter=50, tol=1e-4):
    """Iterate m <- m + gain*(obs - smooth(m)); confine to mask if given."""
    m = obs.copy()
    if mask is not None:
        m = m * mask
    hist = []
    for k in range(1, max_iter + 1):
        r = obs - smooth(m)
        m = m + gain * r
        if mask is not None:
            m = m * mask
        rel = np.max(np.abs(r)) / max(np.max(np.abs(obs)), 1e-30)
        hist.append(rel)
        if rel < tol:
            break
    return m, k, hist


def check_forward():
    truth = cap(15.0, 300.0, 6.0, 1.0)           # 6-deg disc, 1 m EWH
    obs = smooth(truth)
    peak_lost = 1.0 - obs.max()
    print("\nforward modelling")
    print("  filtered peak                : %.4f of truth (%.1f%% lost)"
          % (obs.max(), 100 * peak_lost))

    m, it, hist = forward_model(obs, gain=1.0, max_iter=60, tol=1e-4)
    print("  unconstrained: %2d iters, peak %.4f, residual %.2e"
          % (it, m.max(), hist[-1]))

    mask = (truth > 0).astype(float)
    mm, itm, histm = forward_model(obs, mask=mask, gain=1.0,
                                   max_iter=60, tol=1e-4)
    inside = mask > 0
    rec = mm[inside].mean()
    print("  mask-constrained: %2d iters, mean inside %.4f (truth 1.0)"
          % (itm, rec))
    print("  leaked signal outside the cap: raw %.4f -> corrected %.4f"
          % (np.abs(obs[~inside]).max(), np.abs(mm[~inside]).max()))
    # the corrected field must beat the filtered field inside the cap
    err_raw = np.abs(obs[inside] - 1.0).mean()
    err_cor = np.abs(mm[inside] - 1.0).mean()
    print("  mean |error| inside: raw %.4f -> corrected %.4f (%.1fx better)"
          % (err_raw, err_cor, err_raw / err_cor))
    assert err_cor < err_raw / 3.0
    # monotone convergence of the residual
    assert all(np.diff(histm) < 1e-12), "residual must decrease"
    # a zero field is a fixed point (no spurious mass invented)
    z, _, _ = forward_model(np.zeros_like(obs), max_iter=5)
    assert np.max(np.abs(z)) < 1e-15, "zero must stay zero"
    print("  zero input stays zero        : ok")
    return dict(peak=obs.max(), iters=itm, rec=rec,
                err_raw=err_raw, err_cor=err_cor)


def check_gain():
    """gain > 1 accelerates; too large diverges - the caller needs a bound."""
    truth = cap(15.0, 300.0, 6.0, 1.0)
    obs = smooth(truth)
    mask = (truth > 0).astype(float)
    print("\n  gain sweep (mask-constrained, tol 1e-4, max 60):")
    out = {}
    for g in [0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 9.0]:
        m, it, hist = forward_model(obs, mask=mask, gain=g,
                                    max_iter=60, tol=1e-4)
        div = not np.isfinite(m).all() or hist[-1] > hist[0] \
            or np.max(np.abs(m)) > 1e3
        print("    gain %.1f: %2d iters, residual %.2e%s"
              % (g, it, hist[-1], "  DIVERGING" if div else ""))
        out[g] = (it, hist[-1], div)
    return out


# ----------------------------------------------------- gridded scaling
def grid_scaling(model_series, min_signal=1e-3):
    """Per-pixel least-squares gain k = sum(true*filt)/sum(filt^2).

    Where the model carries (almost) no signal the denominator vanishes
    and k is meaningless - a ratio of two numerical zeros. Those pixels
    must come back NaN, not a large number that silently multiplies the
    data. MIN_SIGNAL is relative to the strongest pixel.
    """
    num = np.zeros_like(model_series[0])
    den = np.zeros_like(model_series[0])
    for tru in model_series:
        f = smooth(tru)
        num += tru * f
        den += f * f
    weak = den < min_signal * den.max()
    k = np.full_like(den, np.nan)
    k[~weak] = num[~weak] / den[~weak]
    return k


def check_scaling():
    rng = np.random.default_rng(3)
    base = cap(15.0, 300.0, 6.0, 1.0)
    side = cap(24.0, 308.0, 5.0, 1.0)          # a neighbouring basin
    inside = base > 0
    amp = 1.0 + 0.4 * np.sin(2 * np.pi * np.arange(24) / 12)
    # the MODEL: the target basin breathing, with an independent
    # neighbour that also leaks into it
    series = [base * a + side * b for a, b in
              zip(amp, 0.6 * amp[::-1] + 0.05 * rng.normal(size=24))]
    k = grid_scaling(series)
    print("\ngridded scaling")
    print("  k inside the target basin : min %.3f  median %.3f  max %.3f"
          % (np.nanmin(k[inside]), np.nanmedian(k[inside]),
             np.nanmax(k[inside])))
    assert np.nanmedian(k[inside]) > 1.0, "a smoother needs k > 1 inside"

    # apply k to a field with a DIFFERENT mixture than the model had -
    # this is the real use and the real caveat, so the improvement must
    # be large but NOT exact
    test = base * 1.7 + side * 0.2
    filt = smooth(test)
    corr = k * filt
    err_raw = np.abs(filt[inside] - 1.7).mean()
    err_cor = np.abs(corr[inside] - 1.7).mean()
    print("  mismatched epoch: mean |error| %.4f -> %.4f (%.1fx better)"
          % (err_raw, err_cor, err_raw / err_cor))
    assert err_cor < err_raw / 2.0, "k must help"
    assert err_cor > 1e-6, ("k reproduced the truth exactly - the test "
                            "pattern is not independent of the model")

    # scale invariance: k depends on the model PATTERN, not its amplitude
    k2 = grid_scaling([f * 1000.0 for f in series])
    d = np.nanmax(np.abs(k2[inside] - k[inside]))
    print("  k is amplitude-invariant  : maxdiff %.2e" % d)
    assert d < 1e-9

    # and the caveat, quantified: a model with the WRONG pattern gives a
    # different k. This is why INFO must report the model provenance.
    kWrong = grid_scaling([side * a for a in amp])
    covered = np.isfinite(kWrong[inside])
    print("  model without signal in the basin: %.0f%% of the target "
          "pixels come back NaN" % (100 * (1 - covered.mean())))
    assert not covered.all(), ("pixels the model does not cover must be "
                              "NaN, not a ratio of numerical zeros")
    if covered.any():
        dk = np.nanmedian(np.abs(kWrong[inside][covered]
                                 - k[inside][covered]))
        print("  where it does overlap, k shifts by %.3f" % dk)
        assert dk > 0.05, "the pattern dependence must be visible"
    return dict(kmed=float(np.nanmedian(k[inside])),
                err_raw=err_raw, err_cor=err_cor)


def check_pinned():
    """Small, exactly reproducible values to pin in the MATLAB suite."""
    print("\npinned reference values (nmax=%d, Gauss 500 km)" % NMAX)
    truth = cap(0.0, 0.0, 10.0, 1.0)
    obs = smooth(truth)
    mask = (truth > 0).astype(float)
    m, it, hist = forward_model(obs, mask=mask, gain=1.0, max_iter=40,
                                tol=1e-6)
    print("  cap(0,0,10deg,1.0): filtered peak      %.10f" % obs.max())
    print("                      iterations to 1e-6 %d" % it)
    print("                      mean inside after  %.10f"
          % m[mask > 0].mean())
    print("                      max outside after  %.10f"
          % np.abs(m[mask == 0]).max())
    print("  gaussian W(1), W(10), W(40)          : %.10f %.10f %.10f"
          % (GW[1], GW[10], GW[40]))


if __name__ == "__main__":
    check_transforms()
    fwd = check_forward()
    check_gain()
    scl = check_scaling()
    check_pinned()
    print("\nvalidate_leakage: all identities hold")
