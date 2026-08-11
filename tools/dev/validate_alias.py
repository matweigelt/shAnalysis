"""validate_alias.py - S2 tidal alias (161 d) removal, Python first.

GravIS fits a harmonic at the S2 alias period together with bias,
trend, annual and semi-annual, and subtracts ONLY that harmonic from
each monthly product (https://gravis.gfz.de/corrections).

The subtlety: across the GRACE / GRACE-FO boundary the nodal planes are
not aligned, so a single harmonic fitted over the joined record is
wrong. Landerer et al. (2020) prescribe a 100 degree PHASE OFFSET
between the two missions. One amplitude pair is still estimated - the
offset is fixed, not free - so the model gains no parameters but does
gain the ability to describe both missions at once.

Checks: exact recovery of a synthetic alias, that ignoring the offset
biases the recovered amplitude, and that only the alias is removed.

Developed by Matthias Weigelt with the help of Claude (Opus 5),
2026-08-11 (v3.4.0).
"""
import numpy as np

P_S2 = 161.0 / 365.25          # years
PHI = np.deg2rad(100.0)


def design(t, split, period=P_S2, phi=PHI):
    """[bias trend cosA sinA cos2A sin2A cosAlias sinAlias]."""
    w = 2 * np.pi
    a = w * t / period + phi * (t >= split)     # the fixed offset
    return np.column_stack([
        np.ones_like(t), t - t.mean(),
        np.cos(w * t), np.sin(w * t),
        np.cos(2 * w * t), np.sin(2 * w * t),
        np.cos(a), np.sin(a)])


def main():
    rng = np.random.default_rng(4)
    # a GRACE/GRACE-FO-like sampling with the 2017-2018 gap
    t = np.concatenate([np.arange(2002.3, 2017.5, 1 / 12),
                        np.arange(2018.5, 2026.0, 1 / 12)])
    split = 2018.0
    truthAmp = np.array([0.7, -0.4])
    y = (2.0 - 1.3 * (t - t.mean())
         + 0.9 * np.cos(2 * np.pi * t) + 0.2 * np.sin(2 * np.pi * t))
    A = design(t, split)
    y = y + A[:, 6] * truthAmp[0] + A[:, 7] * truthAmp[1]
    y = y + 0.01 * rng.normal(size=t.size)

    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    print("with the 100 deg offset : alias amplitude %.4f %.4f (truth %.4f %.4f)"
          % (c[6], c[7], *truthAmp))
    assert np.allclose(c[6:8], truthAmp, atol=0.02), c[6:8]

    # ignoring the offset must degrade the estimate - if it does not, the
    # offset is decorative and the test is not testing anything
    A0 = design(t, split, phi=0.0)
    c0, *_ = np.linalg.lstsq(A0, y, rcond=None)
    err = np.hypot(*(c0[6:8] - truthAmp))
    print("without the offset      : alias amplitude %.4f %.4f (error %.3f)"
          % (c0[6], c0[7], err))
    assert err > 0.1, "the phase offset must matter, or why prescribe it"

    # only the alias is removed: the trend must survive untouched
    resid = y - A[:, 6:8] @ c[6:8]
    A2 = design(t, split)
    c2, *_ = np.linalg.lstsq(A2, resid, rcond=None)
    print("after removal           : trend %.4f (was %.4f), alias %.2e %.2e"
          % (c2[1], c[1], c2[6], c2[7]))
    assert abs(c2[1] - c[1]) < 1e-10, "removing the alias must not touch the trend"
    assert np.allclose(c2[6:8], 0, atol=1e-10), "the alias must be gone"

    # a series entirely on one side of the split is unaffected by the offset
    tg = np.arange(2002.3, 2017.5, 1 / 12)
    Ag = design(tg, split)
    Ag0 = design(tg, split, phi=0.0)
    print("single-mission span     : offset changes design by %.2e"
          % np.max(np.abs(Ag - Ag0)))
    assert np.allclose(Ag, Ag0)
    print("\nvalidate_alias: recovery exact, offset demonstrably matters")


if __name__ == "__main__":
    main()
