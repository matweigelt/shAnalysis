"""Pole-tide convention conversion - numeric pre-validation (v2.7.0).

Mean-pole models (IERS Conventions):
  IERS2010 (Table 7.25), t < 2010.0 cubic, t >= 2010.0 linear [mas]:
    cubic : xm = 55.974 + 1.8243 dt + 0.18413 dt^2 + 0.007024 dt^3
            ym = 346.346 + 1.7896 dt - 0.10729 dt^2 - 0.000908 dt^3
    linear: xm = 23.513 + 7.6141 dt ; ym = 358.891 - 0.6287 dt
  IERS2018 secular pole (2018 update, used by RL06):
    xm = 55.0 + 1.677 dt ; ym = 320.5 + 3.460 dt          [mas]
  dt = t - 2000.0

Wobble params: m1 = xp - xm ; m2 = -(yp - ym)  [arcsec]
Solid pole tide (IERS): dC21 = -1.333e-9 (m1 + 0.0115 m2)
                        dS21 = -1.333e-9 (m2 - 0.0115 m1)
Ocean pole tide (Desai): dC21 = -2.1778e-10 (m1 - 0.01724 m2)
                         dS21 = -2.1778e-10 (m2 + 0.03365 m1)

Converting a solution processed with mean pole A to convention B (same
observed pole xp, yp): the center REMOVED the tide computed with m(A);
under B the removed part should be m(B). Adjustment to the published
coefficients:  dX(A->B) = corr(m(A)) - corr(m(B)) applied with the
formulas above on  dm1 = m1A - m1B = xmB - xmA,
                   dm2 = m2A - m2B = -(ymB - ymA) = ymA - ymB.
"""
import numpy as np

def meanpole(model, t):
    dt = t - 2000.0
    if model == "IERS2010":
        if t < 2010.0:
            xm = 55.974 + 1.8243*dt + 0.18413*dt**2 + 0.007024*dt**3
            ym = 346.346 + 1.7896*dt - 0.10729*dt**2 - 0.000908*dt**3
        else:
            xm = 23.513 + 7.6141*dt
            ym = 358.891 - 0.6287*dt
    elif model == "IERS2018":
        xm = 55.0 + 1.677*dt
        ym = 320.5 + 3.460*dt
    return xm * 1e-3, ym * 1e-3          # mas -> arcsec

def convert(t, A, B, mode="both"):
    xmA, ymA = meanpole(A, t); xmB, ymB = meanpole(B, t)
    dm1 = xmB - xmA                       # m1A - m1B = xmB - xmA
    dm2 = ymA - ymB                       # m2A - m2B = ymA - ymB
    dC = dS = 0.0
    if mode in ("solid", "both"):
        dC += -1.333e-9 * (dm1 + 0.0115*dm2)
        dS += -1.333e-9 * (dm2 - 0.0115*dm1)
    if mode in ("ocean", "both"):
        dC += -2.1778e-10 * (dm1 - 0.01724*dm2)
        dS += -2.1778e-10 * (dm2 + 0.03365*dm1)
    return dC, dS

# identities
for t in [2005.3, 2015.7]:
    dC, dS = convert(t, "IERS2010", "IERS2010")
    assert dC == 0 and dS == 0
    dC1, dS1 = convert(t, "IERS2010", "IERS2018")
    dC2, dS2 = convert(t, "IERS2018", "IERS2010")
    assert abs(dC1 + dC2) < 1e-25 and abs(dS1 + dS2) < 1e-25   # A->B = -(B->A)
print("identities: A->A = 0, A->B = -(B->A)  OK")

# magnitudes + the science: the convention difference is mostly a TREND
for t in [2005.0, 2010.0, 2015.0, 2020.0]:
    dC, dS = convert(t, "IERS2010", "IERS2018")
    print(f"t={t}: dC21={dC:+.3e}  dS21={dS:+.3e}")
r1 = np.array(convert(2020.0, "IERS2010", "IERS2018"))
r0 = np.array(convert(2010.0, "IERS2010", "IERS2018"))
rate = (r1 - r0) / 10
print(f"trend 2010-2020: dC21/yr={rate[0]:+.3e}, dS21/yr={rate[1]:+.3e}")
# published scale: the 2010-linear vs 2018-secular ym rates differ by
# 3.460-(-0.6287)=4.089 mas/yr -> dS21 trend ~ -1.55e-9*(-4.089e-3)
expect_dS_rate = -(1.333e-9 + 2.1778e-10) * (-(3.460 - (-0.6287)) * 1e-3) * -1
# careful: dm2/dt = d(ymA-ymB)/dt = (-0.6287 - 3.460)e-3 = -4.089e-3 "/yr
dm2_rate = (-0.6287 - 3.460) * 1e-3
approx = -(1.333e-9)*(dm2_rate*0 + 0)  # cross terms tiny; dominant:
dS_rate_expect = (-1.333e-9 - 2.1778e-10) * dm2_rate
print(f"expected dS21/yr (dominant term): {dS_rate_expect:+.3e}")
assert abs(rate[1] - dS_rate_expect) / abs(dS_rate_expect) < 0.05
# EWH-relevance sanity: |dC21| ~ 1e-10 corresponds to mm-level C21 EWH - the
# size of published inter-convention C21/S21 trend discrepancies. OK if within
assert 1e-11 < abs(r1[1]) < 1e-9
print("magnitude + dominant-term validation OK")
