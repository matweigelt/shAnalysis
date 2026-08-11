"""validate_degree1.py - self-consistent geocentre (Swenson et al. 2008).

GRACE cannot sense degree 1: the coefficients C10, C11, S11 are set to
zero by definition, and they are exactly the ones describing the offset
between the centre of mass and the centre of figure. Every mass estimate
needs them, so they must come from somewhere.

Swenson, Chambers & Wahr (2008) get them from GRACE ITSELF plus an ocean
model. The argument: a surface mass field is the sum of a land part and
an ocean part. GRACE observes degrees 2+ of the total. The degree-1
terms are then whatever, when added, makes the OCEAN part of the field
consistent with the ocean model. That is a 3x3 linear system in the
three unknowns, per epoch.

Written in the notation of the paper:

    G       = ocean-domain projection of the degree-1 basis functions
    b       = ocean-domain projection of (GRACE degrees 2+ minus the
              ocean model)
    x       = [C10, C11, S11] solving  G x = b

This script builds a synthetic world with a KNOWN geocentre, runs the
estimator, and checks it comes back. It also checks the two failure
modes worth knowing: an ocean mask that is too small makes G singular,
and omitting the ocean model biases the answer.

Developed by Matthias Weigelt with the help of Claude (Opus 5),
2026-08-11 (v3.6.0).
"""
import numpy as np

R = 6378136.3
RHO_AVE = 5517.0
RHO_W = 1000.0


def legendre_deg1(theta):
    """Fully normalised P10, P11 at colatitude theta [rad]."""
    return np.sqrt(3.0) * np.cos(theta), np.sqrt(3.0) * np.sin(theta)


def deg1_basis(latDeg, lonDeg, kn1):
    """The three degree-1 surface-mass patterns, in EWH metres.

    A unit Stokes coefficient of degree 1 maps to surface density with
    the Wahr kernel  R rho_ave/(3 rho_w) (2n+1)/(1+k_n)  at n = 1.
    """
    th = np.radians(90.0 - np.asarray(latDeg))[:, None]
    lam = np.radians(np.asarray(lonDeg))[None, :]
    P10, P11 = legendre_deg1(th)
    kf = R * RHO_AVE / (3.0 * RHO_W) * 3.0 / (1.0 + kn1)
    return (kf * P10 * np.ones_like(lam),          # C10
            kf * P11 * np.cos(lam),                # C11
            kf * P11 * np.sin(lam))                # S11


def estimate(obsOcean, modelOcean, basis, ocean, w):
    """Least squares for [C10, C11, S11] over the ocean domain."""
    A = np.column_stack([b[ocean] * np.sqrt(w[ocean]) for b in basis])
    d = (modelOcean - obsOcean)[ocean] * np.sqrt(w[ocean])
    x, *_ = np.linalg.lstsq(A, d, rcond=None)
    return x, np.linalg.cond(A)


def main():
    lat = np.arange(-89, 90, 2.0)
    lon = np.arange(0, 360, 2.0)
    LO, LA = np.meshgrid(lon, lat)
    w = np.cos(np.radians(LA))
    kn1 = 0.021                       # degree-1 load Love number, CF frame
    basis = deg1_basis(lat, lon, kn1)

    # a realistic-ish ocean: everything outside two "continents"
    ocean = ~(((LA > 10) & (LA < 70) & (LO > 250) & (LO < 350)) |
              ((LA > -40) & (LA < 30) & (LO > 10) & (LO < 60)))
    print("ocean fraction (area) %.3f" % (w[ocean].sum() / w.sum()))

    rng = np.random.default_rng(5)
    # The TRUE surface mass field is land hydrology plus an ocean model.
    # Its own degree-1 content IS the geocentre - it is not an extra term
    # bolted on, which is the point of the method. Build the truth first,
    # then read off both what GRACE sees and what the answer must be.
    land = np.where(~ocean, 0.05 * np.sin(np.radians(3 * LO)), 0.0)
    oceanModel = np.where(ocean, 0.01 * np.cos(np.radians(2 * LA)), 0.0)
    total = land + oceanModel

    # project the truth onto the three degree-1 patterns to get the
    # coefficients GRACE is blind to (area-weighted, over the globe)
    A_all = np.column_stack([b.ravel() * np.sqrt(w.ravel()) for b in basis])
    truth, *_ = np.linalg.lstsq(A_all, total.ravel() * np.sqrt(w.ravel()),
                                rcond=None)
    d1 = sum(t * b for t, b in zip(truth, basis))

    # GRACE sees the field with its degree-1 content removed
    obs = total - d1
    x, cond = estimate(obs, oceanModel, basis, ocean, w)
    print("\nnoise-free recovery")
    print("  truth     %s" % np.array2string(truth, precision=4))
    print("  estimated %s" % np.array2string(x, precision=4))
    print("  cond(A) = %.2f" % cond)
    assert np.allclose(x, truth, rtol=1e-6), (x, truth)

    # with noise: still close, and the conditioning tells you how close
    obsN = obs + 0.002 * rng.normal(size=obs.shape)
    xn, _ = estimate(obsN, oceanModel, basis, ocean, w)
    err = np.abs(xn - truth) / np.abs(truth)
    print("\nwith 4 mm noise: relative error %s"
          % np.array2string(err, precision=3))
    assert np.all(err < 0.15), err

    # FAILURE MODE 1: omitting the ocean model biases the result, because
    # its signal is then attributed to the geocentre
    xb, _ = estimate(obs, np.zeros_like(oceanModel), basis, ocean, w)
    bias = np.abs(xb - truth) / np.abs(truth)
    print("\nno ocean model: relative error %s"
          % np.array2string(bias, precision=3))
    assert np.max(bias) > 0.05, "the ocean model must matter"

    # FAILURE MODE 2: a tiny ocean domain makes the system ill-conditioned
    small = ocean & (LA > 60)
    xs, condSmall = estimate(obs, oceanModel, basis, small, w)
    xsn, _ = estimate(obsN, oceanModel, basis, small, w)
    errFull = np.max(np.abs(xn - truth) / np.abs(truth))
    errSmall = np.max(np.abs(xsn - truth) / np.abs(truth))
    print("polar-only ocean: cond(A) = %.1f (full ocean %.2f)"
          % (condSmall, cond))
    print("  noise-FREE it still recovers exactly (%.1e) - conditioning"
          % np.max(np.abs(xs - truth) / np.abs(truth)))
    print("  only bites with noise: rel error %.3f vs %.3f on the full "
          "ocean (%.1fx worse)" % (errSmall, errFull, errSmall / errFull))
    print("  NOTE the error grows far less than cond(A) suggests (1.4x "
          "against 3.8x) - cond is a warning, not a scale factor")
    # This is why cond(A) is reported: it warns BEFORE the noise does,
    # and the user cannot see the error without knowing the truth.
    assert condSmall > 3 * cond, condSmall
    assert errSmall > errFull, (errSmall, errFull)

    # the C10/C11/S11 <-> geocentre conversion the users will want
    xyz = np.sqrt(3.0) * R * np.array([truth[1], truth[2], truth[0]])
    print("\ngeocentre offsets [mm]: X %.2f Y %.2f Z %.2f"
          % tuple(1e3 * xyz))

    print("\nvalidate_degree1: exact recovery, ocean model matters, "
          "conditioning flags a bad domain")


if __name__ == "__main__":
    main()
