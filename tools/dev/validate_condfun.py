"""Python validation for shLowLevel.buildCondFun (Kvas 2019, Sec. 2.4).

Validates the covariance-conditioning pipeline
    Sigma_tilde = F * ( Z .* T .* (G Sigma G') ) * F'
with G spectral->EWH-grid, F its exact left inverse (F G = I, guaranteed
by the toolbox synthesisMatrix contract A Y = I), Z the region-block
indicator (Kvas eq. 2.117) and T = exp(-psi/psi0) the distance taper
(eq. 2.120).  The Legendre/quadrature exactness itself is asserted in
the MATLAB suite; here the LINEAR-ALGEBRA structure is validated with an
abstract exact quadrature pair.

Checks:
  C1  identity: Z=1, psi0=inf  =>  Sigma_tilde == Sigma exactly
  C2  Z is PSD (Gram matrix of region indicators); T is PD on the
      sphere; Hadamard chain keeps the conditioned matrix PSD
  C3  a singular empirical Sigma (T < P samples) becomes strictly
      positive definite after tapering: rcond improves (Kvas Fig. 2.5)
  C4  region masking really zeroes cross-region covariance in the
      spatial domain

Prepared by Claude (Fable 5), 2026-08-17.
"""
import numpy as np


def make_quadrature_pair(M, P, rng):
    """Abstract exact synthesis/analysis pair: F @ G = I_P, M >= P."""
    G = rng.standard_normal((M, P))
    F = np.linalg.pinv(G)                      # exact left inverse
    assert np.allclose(F @ G, np.eye(P), atol=1e-10)
    return G, F


def sphere_points(M, rng):
    p = rng.standard_normal((3, M))
    p /= np.linalg.norm(p, axis=0)
    psi = np.arccos(np.clip(p.T @ p, -1, 1))
    return p, psi


def condition(Sigma, G, F, Z, T):
    S = G @ Sigma @ G.T
    S = Z * T * S
    out = F @ S @ F.T
    return 0.5 * (out + out.T)


def run():
    rng = np.random.default_rng(1)
    P, M = 30, 80
    G, F = make_quadrature_pair(M, P, rng)
    pts, psi = sphere_points(M, rng)

    # a full-rank PSD spectral covariance
    A = rng.standard_normal((P, 2 * P))
    Sigma = A @ A.T / (2 * P)

    # C1 identity
    Z1 = np.ones((M, M)); T1 = np.ones((M, M))
    out = condition(Sigma, G, F, Z1, T1)
    assert np.max(np.abs(out - Sigma)) < 1e-8 * np.max(np.abs(Sigma)), "C1"
    print("C1  psi0=inf, no regions -> identity            PASS")

    # C2 PSD chain
    region = (pts[2] >= 0).astype(int)         # hemisphere split
    Z = (region[:, None] == region[None, :]).astype(float)
    ev = np.linalg.eigvalsh(Z)
    assert ev.min() > -1e-10, "C2 Z PSD"
    T = np.exp(-psi / 0.6)
    assert np.linalg.eigvalsh(T).min() > 0, "C2 T PD"
    out = condition(Sigma, G, F, Z, T)
    evo = np.linalg.eigvalsh(out)
    assert evo.min() > -1e-8 * evo.max(), "C2 out PSD"
    print("C2  Z PSD, T PD, conditioned matrix PSD          PASS")

    # C3 singular empirical covariance -> PD after taper
    Tsamp = P // 2
    X = rng.standard_normal((P, Tsamp))
    Sig_sing = X @ X.T / Tsamp                 # rank Tsamp < P
    r_raw = np.linalg.eigvalsh(Sig_sing)
    rc_raw = max(r_raw.min(), 0) / r_raw.max()
    out = condition(Sig_sing, G, F, np.ones((M, M)), T)
    r_c = np.linalg.eigvalsh(out)
    rc_c = r_c.min() / r_c.max()
    assert rc_raw < 1e-12 and rc_c > 1e-10, f"C3 {rc_raw:.1e} -> {rc_c:.1e}"
    print(f"C3  rcond {rc_raw:.1e} -> {rc_c:.1e} after taper       PASS")

    # C4 cross-region zeroing in the spatial domain
    S_spat = Z * (G @ Sigma @ G.T)
    i0 = np.where(region == 0)[0][0]
    i1 = np.where(region == 1)[0][0]
    assert S_spat[i0, i1] == 0.0, "C4"
    print("C4  cross-region spatial covariance zeroed       PASS")

    print("\nAll 4 checks PASS")


if __name__ == "__main__":
    run()
