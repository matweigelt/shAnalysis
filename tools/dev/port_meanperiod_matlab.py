#!/usr/bin/env python3
"""Port of the new shSeries.mean body; checks every MATLAB test expectation.

Developed by Matthias Weigelt with the help of Claude (Opus 5), 2026-08-20.

The MATLAB bridge is unavailable this session, so the algorithm is
re-implemented here exactly as written in shSeries.m - same flatten order
(column-major), same design matrix, same integral weights - and every
number asserted by tests/testCorrectness.m is reproduced.
"""
import numpy as np

FAIL = []


def check(name, got, exp, tol):
    g, e = np.asarray(got), np.asarray(exp)
    ok = np.all(g == e) if g.dtype == bool else np.all(np.abs(g - e) <= tol)
    print(f"  {'ok ' if ok else 'FAIL'} {name:38s} got {np.asarray(got).ravel()[:3]}"
          f"  exp {np.asarray(exp).ravel()[:3]}")
    if not ok:
        FAIL.append(name)


def signal(t):
    return (1.0 + 0.30 * (t - 2004) + 1.20 * np.cos(2 * np.pi * (t - 0.08))
            + 0.35 * np.cos(4 * np.pi * (t - 0.55)))


def sincpi(x):
    x = np.asarray(x, float)
    y = np.ones_like(x)
    nz = x != 0
    y[nz] = np.sin(np.pi * x[nz]) / (np.pi * x[nz])
    return y


def truth(t1, t2):
    tm, d = 0.5 * (t1 + t2), t2 - t1
    return (1.0 + 0.30 * (tm - 2004)
            + 1.20 * np.cos(2 * np.pi * (tm - 0.08)) * sincpi(d / 1.0)
            + 0.35 * np.cos(4 * np.pi * (tm - 0.55)) * sincpi(d / 0.5))


def months(a, b):
    return a + (np.arange(round((b - a) * 12)) + 0.5) / 12.0


M = np.array([[1., 0, 0], [2, 3, 0], [4, 5, 6]])
MS = np.array([[0., 0, 0], [0, 1, 0], [0, 2, 3]])


def series(t):
    """Cs, Ss as (3,3,T) - MATLAB layout."""
    y = signal(t)
    return M[:, :, None] * y[None, None, :], MS[:, :, None] * y[None, None, :]


def flatten(Cs, Ss):
    """MATLAB: [reshape(Cs,Nc,T); reshape(Ss,Nc,T)] - column-major."""
    n1, _, T = Cs.shape
    Nc = n1 * n1
    return np.vstack([Cs.reshape(Nc, T, order="F"), Ss.reshape(Nc, T, order="F")])


def unflat(x, n1):
    Nc = n1 * n1
    return (x[:Nc].reshape(n1, n1, order="F"), x[Nc:].reshape(n1, n1, order="F"))


def mean_matlab(Cs, Ss, ep, rng=None, estimator="arithmetic", periods=()):
    """Mirror of the MATLAB method, arithmetic and model branches."""
    ep = np.asarray(ep, float)
    if rng is not None:
        keep = np.flatnonzero((ep >= rng[0]) & (ep <= rng[1]))
        assert keep.size, "emptyRange"
        Cs, Ss, ep = Cs[:, :, keep], Ss[:, :, keep], ep[keep]
    else:
        rng = [np.nanmin(ep), np.nanmax(ep)]
    atten = np.zeros((0, 2))
    if estimator == "arithmetic":
        Cm, Sm = np.nanmean(Cs, 2), np.nanmean(Ss, 2)
        nC = np.sum(~np.isnan(Cs), 2)
        nS = np.sum(~np.isnan(Ss), 2)
        with np.errstate(invalid="ignore"):
            sc = np.nanstd(Cs, 2, ddof=1) / np.sqrt(np.maximum(nC, 1))
            ss = np.nanstd(Ss, 2, ddof=1) / np.sqrt(np.maximum(nS, 1))
        sc = np.where(nC < 2, np.nan, sc)
        ss = np.where(nS < 2, np.nan, ss)
        epoch = np.nanmean(ep)
    else:
        keep = ~np.any(~np.isfinite(Cs), (0, 1)) & ~np.any(~np.isfinite(Ss), (0, 1))
        Cs, Ss, ep = Cs[:, :, keep], Ss[:, :, keep], ep[keep]
        nPar = 6 + 2 * len(periods)
        tm, D = 0.5 * (rng[0] + rng[1]), rng[1] - rng[0]
        X = flatten(Cs, Ss)                       # P x T
        tc = ep - tm
        cols = [np.ones_like(tc), tc, np.cos(2 * np.pi * tc), np.sin(2 * np.pi * tc),
                np.cos(4 * np.pi * tc), np.sin(4 * np.pi * tc)]
        for p in periods:
            cols += [np.cos(2 * np.pi * tc / p), np.sin(2 * np.pi * tc / p)]
        A = np.column_stack(cols)                  # T x ncol
        coef, *_ = np.linalg.lstsq(A, X.T, rcond=None)     # ncol x P
        res = X.T - A @ coef
        rvar = np.sum(res ** 2, 0) / (len(ep) - nPar)
        pers = np.array([1.0, 0.5] + list(periods))
        s = sincpi(D / pers)
        u = np.zeros(nPar)
        u[0] = 1.0
        u[2::2] = s                                # MATLAB u(3:2:end)
        Xm = u @ coef
        q = u @ np.linalg.solve(A.T @ A, u)
        Xs = np.sqrt(np.maximum(rvar, 0) * max(q, 0))
        n1 = Cs.shape[0]
        Cm, Sm = unflat(Xm, n1)
        sc, ss = unflat(Xs, n1)
        epoch = tm
        atten = np.column_stack([pers, s])
    # info
    used = np.sort(ep[np.isfinite(ep)])
    info = dict(nUsed=used.size, epoch=epoch, attenuation=atten,
                span=[used[0], used[-1]],
                centroidOffset=used.mean() - 0.5 * (rng[0] + rng[1]))
    if used.size > 1:
        dt = np.median(np.diff(used))
        info["coverage"] = min(1, used.size * dt / max(rng[1] - rng[0], 2.2e-16))
        info["maxGapYears"] = np.max(np.diff(used))
    return Cm, Sm, sc, ss, info


# ------------------------------------------- testMeanPeriodGappyWindow
print("testMeanPeriodGappyWindow")
t = months(2004, 2010)
t = t[~((t > 2006.4) & (t < 2006.9))]
Cs, Ss = series(t)
tr = truth(2004, 2010)
check("numel(t)", len(t), 66, 0)
check("truth", tr, 1.9, 1e-12)
Cm, Sm, sc, ss, iA = mean_matlab(Cs, Ss, t, [2004, 2010])
check("arith C(1,1)-truth", Cm[0, 0] - tr, 0.069179522165, 1e-9)
check("info.nUsed", iA["nUsed"], 66, 0)
check("info.coverage", iA["coverage"], 0.916666666667, 1e-9)
check("info.maxGapYears", iA["maxGapYears"], 0.583333333333, 1e-9)
check("info.centroidOffset", iA["centroidOffset"], 0.030303030303, 1e-9)
check("arith epoch = centroid", iA["epoch"], t.mean(), 1e-12)
Cm, Sm, sc, ss, iM = mean_matlab(Cs, Ss, t, [2004, 2010], "model")
check("model C(1,1)", Cm[0, 0], tr, 1e-10)
check("model C(3,2) = 5*truth", Cm[2, 1], 5 * tr, 1e-9)
check("model S(3,3) = 3*truth", Sm[2, 2], 3 * tr, 1e-9)
check("model epoch = centre", iM["epoch"], 2007.0, 1e-12)
check("attenuation periods", iM["attenuation"][:, 0], [1, 0.5], 1e-12)
check("attenuation factors", iM["attenuation"][:, 1], [0, 0], 1e-12)
check("model sigma ~ 0", np.max(np.abs(sc)), 0.0, 1e-8)

# --------------------------------------------- testMeanPeriodAttenuation
print("testMeanPeriodAttenuation")
t = months(2004, 2009.5)
Cs, Ss = series(t)
tr = truth(2004, 2009.5)
check("truth", tr, 1.858457517928, 1e-9)
Cm, *_ = mean_matlab(Cs, Ss, t, [2004, 2009.5])
check("arith - truth", Cm[0, 0] - tr, 0.000385268670, 1e-9)
Cm, Sm, sc, ss, iM = mean_matlab(Cs, Ss, t, [2004, 2009.5], "model")
check("model C(1,1)", Cm[0, 0], tr, 1e-10)
check("attenuation D=5.5", iM["attenuation"][:, 1], [-0.057874524761, 0], 1e-9)

# ------------------------------------ testMeanRangeMatchesSelectThenMean
print("testMeanRangeMatchesSelectThenMean")
t = months(2004, 2010)
Cs, Ss = series(t)
a = mean_matlab(Cs, Ss, t, [2005, 2008])[0]
k = (t >= 2005) & (t <= 2008)
b = mean_matlab(Cs[:, :, k], Ss[:, :, k], t[k])[0]
check("Range == select+mean", a, b, 1e-14)
check("full-series mean unchanged", mean_matlab(Cs, Ss, t)[0], Cs.mean(2), 1e-14)

# ----------------------------------- testMeanSigmaUsesValidEpochCount
print("testMeanSigmaUsesValidEpochCount")
t = months(2004, 2007)
Cs, Ss = series(t)
Cs[:, :, [3, 7, 11]] = np.nan
Ss[:, :, [3, 7, 11]] = np.nan
Cm, Sm, sc, ss, _ = mean_matlab(Cs, Ss, t)
check("C(1,1)", Cm[0, 0], 1.473863636364, 1e-9)
check("sigmaC(1,1)", sc[0, 0], 0.156294825874, 1e-9)
v = Cs[0, 0, :]
v = v[~np.isnan(v)]
check("valid count", v.size, 33, 0)
check("sigma = std/sqrt(33)", sc[0, 0], v.std(ddof=1) / np.sqrt(33), 1e-12)
Cs2 = Cs.copy()
Cs2[0, 0, :] = np.nan
Cs2[0, 0, 0] = 5.0
Cm2, _, sc2, _, _ = mean_matlab(Cs2, Ss, t)
check("single-sample mean", Cm2[0, 0], 5.0, 1e-12)
check("single-sample sigma is NaN", np.isnan(sc2[0, 0]), True, 0)

print("\n%d failure(s)" % len(FAIL), FAIL)
