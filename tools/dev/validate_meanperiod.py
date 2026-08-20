#!/usr/bin/env python3
"""validate_meanperiod: reference for shSeries.mean(Range=, Estimator=).

Python-first validation of the two period-mean estimators before the
MATLAB implementation, and the source of the frozen numbers checked by
testCorrectness/testMeanPeriod*.

The question: given a spherical harmonic series and a user window
[t1, t2], what is the mean field over that window?

  arithmetic  equal weights on the epochs that fall inside the window.
              Unbiased only if the sampling is balanced with respect to
              trend and annual signal inside the window; the estimate
              belongs to the epoch CENTROID, not to the window centre.

  model       least-squares bias/trend/annual/semiannual (+Periods) fit
              on the epochs inside the window with T0 = the window
              centre, then the EXACT time integral of that model over
              [t1, t2].  Because T0 is the window centre the integral
              collapses to

                  mean = bias + sum_k s(P_k) * cos-coefficient(P_k),
                  s(P)  = sin(pi*D/P) / (pi*D/P),   D = t2 - t1,

              i.e. the model evaluated at the window centre with every
              harmonic damped by a sinc factor.  For a whole number of
              years s(1) = s(1/2) = 0 exactly: the period mean is then
              the trend line at the window centre, whatever the phase of
              the annual signal or where the gaps sit.

A third candidate, Voronoi time-weighting, was tested and REJECTED - see
the table below: it never wins and is 4x worse than the arithmetic mean
across the GRACE/GRACE-FO mission gap, because the two months bracketing
the gap each collect half a year of weight.

Stdlib + numpy only.  Nonzero exit if any frozen value drifts.

Developed by Matthias Weigelt with the help of Claude (Opus 5),
2026-08-20.
"""
import sys

import numpy as np

# ----------------------------------------------------------------- model
# one Stokes coefficient, exactly inside the fit space:
#   x(t) = B + TR*(t-2004) + Aa*cos(2pi(t-pa)) + As*cos(4pi(t-ps))
B, TR, AA, PA, AS, PS = 1.0, 0.30, 1.20, 0.08, 0.35, 0.55


def signal(t):
    return (B + TR * (t - 2004.0)
            + AA * np.cos(2 * np.pi * (t - PA))
            + AS * np.cos(4 * np.pi * (t - PS)))


def truth(t1, t2):
    """Exact analytic mean of `signal` over [t1, t2].

    int cos(w(t-p)) dt / D = cos(w(tm-p)) * sin(w*D/2)/(w*D/2)
    """
    tm, d = 0.5 * (t1 + t2), t2 - t1
    out = B + TR * (tm - 2004.0)
    for amp, per, ph in ((AA, 1.0, PA), (AS, 0.5, PS)):
        w = 2 * np.pi / per
        out += amp * np.cos(w * (tm - ph)) * sinc(d / per)
    return out


def sinc(x):
    """sin(pi*x)/(pi*x), = 1 at x = 0 (NOT the Signal Toolbox function)."""
    x = np.asarray(x, float)
    y = np.ones_like(x)
    nz = x != 0
    y[nz] = np.sin(np.pi * x[nz]) / (np.pi * x[nz])
    return y


# ------------------------------------------------------------ estimators
def design(t, t0, periods=()):
    tc = t - t0
    cols = [np.ones_like(tc), tc,
            np.cos(2 * np.pi * tc), np.sin(2 * np.pi * tc),
            np.cos(4 * np.pi * tc), np.sin(4 * np.pi * tc)]
    for p in periods:
        cols += [np.cos(2 * np.pi * tc / p), np.sin(2 * np.pi * tc / p)]
    return np.column_stack(cols)


def integral_weights(d, periods=()):
    """u such that u @ coef is the exact model mean over a window of
    length d, when T0 is the window centre."""
    u = [1.0, 0.0]
    for p in (1.0, 0.5) + tuple(periods):
        u += [float(sinc(d / p)), 0.0]
    return np.array(u)


def mean_model(t, y, t1, t2, periods=()):
    tm, d = 0.5 * (t1 + t2), t2 - t1
    A = design(t, tm, periods)
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    u = integral_weights(d, periods)
    res = y - A @ coef
    rvar = res @ res / (len(t) - A.shape[1])
    q = u @ np.linalg.solve(A.T @ A, u)
    return u @ coef, np.sqrt(rvar * q)


def mean_arith(y):
    return float(np.mean(y))


def voronoi_weights(t, t1, t2):        # rejected, kept to document why
    e = np.empty(len(t) + 1)
    e[1:-1] = 0.5 * (t[:-1] + t[1:])
    e[0], e[-1] = t1, t2
    w = np.diff(e)
    return w / w.sum()


# ------------------------------------------------------------- scenarios
def months(a, b):
    return a + (np.arange(round((b - a) * 12)) + 0.5) / 12.0


def scenario(name, t1, t2, cut=None):
    t = months(t1, t2)
    if cut is not None:
        t = t[~((t > cut[0]) & (t < cut[1]))]
    return name, t1, t2, t


SCENARIOS = [
    scenario("complete 2004-2010", 2004.0, 2010.0),
    scenario("6-month block missing", 2004.0, 2010.0, (2006.4, 2006.9)),
    scenario("11-month mission gap", 2015.0, 2021.0, (2017.5, 2018.5)),
    scenario("non-integer window", 2004.0, 2009.5),
]

# frozen: deviation of each estimator from the exact period mean.
# testCorrectness reproduces ARITH and MODEL from MATLAB.
FROZEN = {
    "complete 2004-2010":    (+0.000000, +0.000000, +0.000000),
    "6-month block missing": (+0.069180, +0.074579, +0.000000),
    "11-month mission gap":  (+0.000000, -0.128419, +0.000000),
    "non-integer window":    (+0.000385, +0.000385, +0.000000),
}


def main():
    bad = 0
    print("deviation from the exact period mean  [coefficient units]\n")
    print(f"{'scenario':24s} {'n':>3s} {'truth':>9s} {'arith':>10s}"
          f" {'voronoi':>10s} {'model':>10s}")
    for name, t1, t2, t in SCENARIOS:
        y = signal(t)
        tm = truth(t1, t2)
        got = (mean_arith(y) - tm,
               voronoi_weights(t, t1, t2) @ y - tm,
               mean_model(t, y, t1, t2)[0] - tm)
        print(f"{name:24s} {len(t):3d} {tm:9.4f} "
              + " ".join(f"{g:+10.6f}" for g in got))
        exp = FROZEN[name]
        for k, lbl in enumerate(("arith", "voronoi", "model")):
            if abs(got[k] - exp[k]) > 5e-6:
                print(f"  DRIFT {name}/{lbl}: {got[k]:+.6f} != {exp[k]:+.6f}")
                bad += 1

    # the sinc attenuation actually applied, reported in info.attenuation
    print("\nharmonic attenuation s(P) = sin(pi*D/P)/(pi*D/P):")
    for d in (6.0, 5.5, 1.0, 0.5):
        s = sinc(np.array([d / 1.0, d / 0.5]))
        print(f"  D = {d:4.2f} yr   annual {s[0]:+.6f}   semiannual {s[1]:+.6f}")
    if abs(sinc(np.array([6.0]))[0]) > 1e-12:
        print("  DRIFT: whole-year window must annihilate the annual term")
        bad += 1
    # frozen for testMeanPeriodAttenuation (D = 5.5 yr)
    s55 = sinc(np.array([5.5, 11.0]))
    print(f"  frozen D = 5.5: {s55[0]:.9f} {s55[1]:.9f}")
    for got_, exp_ in zip(s55, (-0.057874525, 0.0)):
        if abs(got_ - exp_) > 1e-8:
            print(f"  DRIFT sinc: {got_:.9f} != {exp_:.9f}")
            bad += 1

    # standard error of the arithmetic mean must divide by the number of
    # VALID epochs, not by nEpochs - the bug this PR fixes
    print("\nstandard error with NaN epochs:")
    y = signal(months(2004.0, 2007.0))
    y[[3, 7, 11]] = np.nan
    n_valid = np.count_nonzero(~np.isnan(y))
    se_right = np.nanstd(y, ddof=1) / np.sqrt(n_valid)
    se_wrong = np.nanstd(y, ddof=1) / np.sqrt(len(y))
    print(f"  n = {len(y)}, valid = {n_valid}: correct {se_right:.6f}, "
          f"old (÷sqrt(nEpochs)) {se_wrong:.6f} "
          f"-> understated by {100 * (1 - se_wrong / se_right):.1f}%")

    print("\nvalidate_meanperiod: %d drift(s)" % bad)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
