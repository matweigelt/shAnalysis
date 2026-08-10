import numpy as np

def legendre_alf(nmax, lat):
    """Fully normalized (4-pi) associated Legendre functions, same recursion
    as legendreALF.m, implemented independently in Python for cross-check."""
    lat = np.atleast_1d(lat)
    t = np.sin(lat)
    u = np.cos(lat)
    nlat = lat.size
    P = np.zeros((nmax+1, nmax+1, nlat))
    P[0,0,:] = 1.0
    if nmax >= 1:
        P[1,0,:] = np.sqrt(3)*t
        P[1,1,:] = np.sqrt(3)*u
    for n in range(2, nmax+1):
        P[n,n,:] = u*np.sqrt((2*n+1)/(2*n))*P[n-1,n-1,:]
        P[n,n-1,:] = t*np.sqrt(2*n+1)*P[n-1,n-1,:]
        for m in range(0, n-1):
            a = np.sqrt(((2*n-1)*(2*n+1))/((n-m)*(n+m)))
            b = np.sqrt(((2*n+1)*(n+m-1)*(n-m-1))/((2*n-3)*(n-m)*(n+m)))
            P[n,m,:] = a*t*P[n-1,m,:] - b*P[n-2,m,:]
    return P

# ---- Test 1: orthonormality on the sphere via Gauss-Legendre quadrature ----
nmax = 40
nodes, weights = np.polynomial.legendre.leggauss(nmax+5)   # nodes = cos(colat) essentially; use as sin(lat)
lat = np.arcsin(nodes)   # sin(lat) = nodes
P = legendre_alf(nmax, lat)

# integral over full sphere of Pbar_nm^2 * (2 - delta_m0) dOmega should be 4*pi,
# and cross-degree/order integrals should vanish.
# dOmega = dlon dlat*cos(lat) -> integrating lon analytically gives 2*pi (or pi for cross m sin/cos),
# and the Gauss-Legendre weights already handle the d(sin lat) = cos(lat) dlat integration.
max_err_diag = 0.0
max_err_offdiag = 0.0
for n in range(nmax+1):
    for m in range(n+1):
        lambda_factor = 2*np.pi if m == 0 else np.pi
        integral = np.sum(weights * P[n,m,:]**2) * lambda_factor
        expected = 4*np.pi
        max_err_diag = max(max_err_diag, abs(integral/expected - 1))

# a few off-diagonal (different degree, same order) orthogonality checks
for (n1, n2, m) in [(2,4,0), (3,5,1), (10,12,3), (20,25,5)]:
    integral = np.sum(weights * P[n1,m,:]*P[n2,m,:]) * 2*np.pi
    max_err_offdiag = max(max_err_offdiag, abs(integral))

print(f"Test 1 - orthonormality: max relative error (diagonal) = {max_err_diag:.3e}")
print(f"Test 1 - orthogonality: max abs error (off-diagonal)   = {max_err_offdiag:.3e}")

# ---- Test 2: round-trip synthesis -> quadrature analysis recovers coefficients ----
np.random.seed(42)
nmax2 = 20
C = np.zeros((nmax2+1, nmax2+1))
S = np.zeros((nmax2+1, nmax2+1))
for n in range(nmax2+1):
    for m in range(n+1):
        C[n,m] = np.random.randn() * 1e-6
        if m > 0:
            S[n,m] = np.random.randn() * 1e-6

nlat_q = nmax2 + 5
nodes2, weights2 = np.polynomial.legendre.leggauss(nlat_q)
lat_q = np.arcsin(nodes2)
nlon_q = 2*(2*nmax2+1)
lon_q = np.linspace(0, 2*np.pi, nlon_q, endpoint=False)

Pq = legendre_alf(nmax2, lat_q)

grid = np.zeros((lat_q.size, lon_q.size))
for k in range(lat_q.size):
    for n in range(nmax2+1):
        Pn = Pq[n, :n+1, k]
        for m in range(n+1):
            grid[k,:] += Pn[m]*(C[n,m]*np.cos(m*lon_q) + S[n,m]*np.sin(m*lon_q))

# quadrature analysis: recover Cnm, Snm
Chat = np.zeros_like(C)
Shat = np.zeros_like(S)
dlon = 2*np.pi/nlon_q
for n in range(nmax2+1):
    for m in range(n+1):
        Pnm = Pq[n,m,:]
        integrand_c = grid * np.cos(m*lon_q)[None,:]
        integrand_s = grid * np.sin(m*lon_q)[None,:]
        lon_int_c = np.sum(integrand_c, axis=1) * dlon
        lon_int_s = np.sum(integrand_s, axis=1) * dlon
        Chat[n,m] = np.sum(weights2 * Pnm * lon_int_c) / (4*np.pi)
        Shat[n,m] = np.sum(weights2 * Pnm * lon_int_s) / (4*np.pi)

errC = np.max(np.abs(Chat - C))
errS = np.max(np.abs(Shat[:, 1:] - S[:, 1:])) if nmax2 > 0 else 0.0
print(f"Test 2 - round-trip max abs error in Cnm: {errC:.3e}")
print(f"Test 2 - round-trip max abs error in Snm: {errS:.3e}")

assert max_err_diag < 1e-8, "Legendre normalization failed"
assert max_err_offdiag < 1e-8, "Legendre orthogonality failed"
assert errC < 1e-15 * 10, "Round-trip synthesis/analysis failed for Cnm"
assert errS < 1e-15 * 10, "Round-trip synthesis/analysis failed for Snm"
print("\nAll validation checks passed.")
