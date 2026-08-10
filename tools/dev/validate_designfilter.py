"""designFilter numerics: W = (N + a*inv(S))^-1 N per order block."""
import numpy as np
# scalar/diagonal gain: sigma=3e-11 noise, Kaula K/n^2 signal RMS, alpha=1
K = 1e-5 * 1e-5  # placeholder scale; formula is scale-covariant in ratio
def gain(sig, s, a): return 1.0 / (1.0 + a * sig**2 / s**2)
sig = 3e-11
for n in [10, 30, 60]:
    s = 1e-6 / n**2
    print(f"n={n}: gain={gain(sig, s, 1.0):.6f}")
# block case: SPD N, SPD S -> eigenvalues of W in (0,1); W S-symmetric
rng = np.random.default_rng(5)
A = rng.normal(size=(4, 4)); N = A @ A.T + 4*np.eye(4)
B = rng.normal(size=(4, 4)); S = B @ B.T + np.eye(4)
a = 0.7
W = np.linalg.solve(N + a*np.linalg.inv(S), N)
ev = np.linalg.eigvals(W)
assert np.all(np.isreal(ev)) or np.max(np.abs(ev.imag)) < 1e-12
ev = ev.real
print("block eigenvalues:", np.round(np.sort(ev), 4))
assert np.all(ev > 0) and np.all(ev < 1)
# limits
W0 = np.linalg.solve(N + 1e-12*np.linalg.inv(S), N)
assert np.allclose(W0, np.eye(4), atol=1e-9)          # alpha->0: identity
# pinned 2x2 case for the MATLAB test
N2 = np.array([[4.0, 1.0], [1.0, 3.0]])
S2 = np.array([[2.0, 0.5], [0.5, 1.0]])
W2 = np.linalg.solve(N2 + 0.5*np.linalg.inv(S2), N2)
np.set_printoptions(precision=15)
print("W2 =", W2.flatten())
print("ALL OK")
