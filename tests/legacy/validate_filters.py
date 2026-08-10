import numpy as np
from scipy.special import eval_legendre
from scipy.integrate import quad

def jekeli_Wn_recursion(nmax, radius_km, R_km=6378.1363):
    b = np.log(2) / (1 - np.cos(radius_km / R_km))
    W = np.zeros(nmax+1)
    W[0] = 1.0
    if nmax >= 1:
        W[1] = (1+np.exp(-2*b))/(1-np.exp(-2*b)) - 1/b
    for n in range(1, nmax):
        W[n+1] = -(2*n+1)/b * W[n] + W[n-1]
    return W, b

nmax = 60
radius_km = 300.0
R_km = 6378.1363
W_rec, b = jekeli_Wn_recursion(nmax, radius_km, R_km)

def Wspatial(x):
    return b/(2*np.pi) * np.exp(-b*(1-x)) / (1 - np.exp(-2*b))

max_err = 0.0
for n in [0, 1, 2, 5, 10, 20, 40, 60]:
    integrand = lambda x: Wspatial(x) * eval_legendre(n, x)
    val, err = quad(integrand, -1, 1, limit=200)
    Wn_quad = 2*np.pi*val
    diff = abs(Wn_quad - W_rec[n])
    max_err = max(max_err, diff)
    print(f"n={n:3d}  recursion={W_rec[n]: .6e}  quadrature={Wn_quad: .6e}  diff={diff:.2e}")

print(f"\nmax abs diff recursion vs. direct spatial-kernel quadrature: {max_err:.3e}")
assert max_err < 1e-6, "Jekeli Gaussian filter recursion mismatch"
print("Gaussian filter (Jekeli 1981) recursion validated against spatial-domain kernel.")
