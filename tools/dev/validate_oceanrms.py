"""validate_oceanrms.py - the open-ocean noise metric, Python first.

GRACE processing centres quote the RMS of a field over the OPEN ocean -
ocean points more than ~1000 km from any coast - as the standard noise
metric (Dahle et al. 2025, Sect. 2.1.1; also the uncertainty basis for
the GravIS TWS and OBP products). It is the natural estimate for
leakageCorrect's NoiseLevel, because far from land a GRACE field should
contain almost no real signal: what is left is error.

Two things have to be right:
  1. the EROSION - "more than d km from any coast" is the ocean mask
     shrunk by d, i.e. every point whose great-circle distance to the
     nearest non-ocean point exceeds d;
  2. the AVERAGING - grid cells are not equal in area, so an unweighted
     RMS over a lat/lon grid counts polar cells far too heavily.

Checks the erosion against an analytic spherical cap and the weighting
against a field whose area-weighted RMS is known.

Developed by Matthias Weigelt with the help of Claude (Opus 5),
2026-08-11 (v3.5.1).
"""
import numpy as np

R = 6371.0


def gc(la1, lo1, la2, lo2):
    """great-circle distance [km], radians in."""
    return R * np.arccos(np.clip(
        np.sin(la1) * np.sin(la2) +
        np.cos(la1) * np.cos(la2) * np.cos(lo1 - lo2), -1, 1))


def erode(ocean, LA, LO, dkm):
    """Ocean points farther than dkm from the nearest non-ocean point."""
    if dkm <= 0:
        return ocean.copy()
    la = np.radians(LA[ocean]); lo = np.radians(LO[ocean])
    nla = np.radians(LA[~ocean]); nlo = np.radians(LO[~ocean])
    if nla.size == 0:
        return ocean.copy()
    dmin = np.full(la.size, np.inf)
    for c in range(0, nla.size, 500):
        s = slice(c, c + 500)
        dmin = np.minimum(dmin, gc(la[:, None], lo[:, None],
                                   nla[s][None, :], nlo[s][None, :]).min(1))
    out = np.zeros_like(ocean)
    out[ocean] = dmin > dkm
    return out


def main():
    lat = np.arange(-89, 90, 2.0)
    lon = np.arange(0, 360, 2.0)
    LO, LA = np.meshgrid(lon, lat)

    # an "ocean" = everything outside a 30 deg cap at the north pole
    psi = np.degrees(np.arccos(np.clip(np.sin(np.radians(90.0)) *
                    np.sin(np.radians(LA)), -1, 1)))
    ocean = psi > 30.0
    print("cap-land model: ocean fraction (area) %.4f"
          % (np.cos(np.radians(LA))[ocean].sum() /
             np.cos(np.radians(LA)).sum()))

    # eroding by d must move the boundary by d: the surviving points
    # should start at colatitude 30 deg + d/R (in degrees)
    for dkm in (0.0, 500.0, 1000.0, 2000.0):
        e = erode(ocean, LA, LO, dkm)
        edge = psi[e].min()
        want = 30.0 + np.degrees(dkm / R)
        print("  erode %5.0f km: boundary at %.2f deg (analytic %.2f, "
              "grid step 2.0)" % (dkm, edge, want))
        assert abs(edge - want) <= 2.5, (edge, want)

    # area weighting: a field that is 1 north of the equator and -1 south
    # has RMS 1 either way, but a field weighted by cos(lat) must not be
    # dominated by the polar rows. Use a field that IS latitude dependent.
    f = np.cos(np.radians(LA))                      # small at the poles
    w = np.cos(np.radians(LA))
    rms_w = np.sqrt((w * f ** 2).sum() / w.sum())
    rms_u = np.sqrt((f ** 2).mean())
    print("\nweighting on a cos(lat) field: weighted %.4f, unweighted %.4f"
          % (rms_w, rms_u))
    assert rms_w > rms_u, "unweighted over-counts the small polar values"
    # the analytic area-weighted value of cos^2 over the sphere is 2/3
    print("  weighted^2 = %.4f (analytic mean of cos^2 = %.4f)"
          % (rms_w ** 2, 2.0 / 3.0))
    assert abs(rms_w ** 2 - 2.0 / 3.0) < 0.01

    # and on white noise both agree, which is why the difference is easy
    # to miss until the field has structure
    rng = np.random.default_rng(1)
    g = rng.normal(size=LA.shape)
    a = np.sqrt((w * g ** 2).sum() / w.sum())
    b = np.sqrt((g ** 2).mean())
    print("  white noise: weighted %.4f, unweighted %.4f (agree)" % (a, b))
    assert abs(a - b) < 0.05

    print("\nvalidate_oceanrms: erosion matches the analytic boundary, "
          "area weighting matches cos^2")


if __name__ == "__main__":
    main()
