"""Numeric pre-validation for the v2.6.0 comparison suite."""
import numpy as np
from scipy import stats
rng = np.random.default_rng(7)

# ---- 1. degree correlation + diff spectrum identities
L = 20; n1 = L + 1
C1 = np.tril(rng.normal(size=(n1, n1))); S1 = np.tril(rng.normal(size=(n1, n1)), -0)
S1[:, 0] = 0
C2 = 0.5 * C1; S2 = 0.5 * S1          # pure scaling
def degcorr(Ca, Sa, Cb, Sb):
    r = np.zeros(n1)
    for n in range(n1):
        num = (Ca[n, :n+1] * Cb[n, :n+1]).sum() + (Sa[n, :n+1] * Sb[n, :n+1]).sum()
        da = (Ca[n, :n+1]**2 + Sa[n, :n+1]**2).sum()
        db = (Cb[n, :n+1]**2 + Sb[n, :n+1]**2).sum()
        r[n] = num / np.sqrt(da * db) if da * db > 0 else np.nan
    return r
r = degcorr(C1, S1, C2, S2)
assert np.allclose(r[2:], 1.0), "scaled field must have degree corr 1"
dAmp = np.sqrt([( (C1-C2)[n,:n+1]**2 + (S1-S2)[n,:n+1]**2 ).sum() for n in range(n1)])
sAmp = np.sqrt([( C1[n,:n+1]**2 + S1[n,:n+1]**2 ).sum() for n in range(n1)])
assert np.allclose(dAmp, 0.5 * sAmp), "diff amp of half-field = half signal amp"
print("1. degree corr + diff spectrum: OK")

# ---- 2. weighted spatial stats + Taylor identity
lat = np.arange(-88, 89, 4.0); lon = np.arange(0, 360, 6.0)
w = np.cos(np.deg2rad(lat))[:, None] * np.ones((1, lon.size)); w /= w.sum()
A = rng.normal(size=(lat.size, lon.size)); B = 0.8 * A + 0.3 * rng.normal(size=A.shape)
wm = lambda X: (w * X).sum()
a = A - wm(A); b = B - wm(B)
stdA = np.sqrt(wm(a*a)); stdB = np.sqrt(wm(b*b))
corr = wm(a*b) / (stdA * stdB)
crmsd = np.sqrt(wm((a-b)**2))
assert abs(crmsd**2 - (stdA**2 + stdB**2 - 2*stdA*stdB*corr)) < 1e-12
print(f"2. Taylor identity: OK (corr={corr:.4f}, ratio={stdB/stdA:.4f})")

# ---- 3. NSE
ref = rng.normal(size=200); mod = ref + 0.5 * rng.normal(size=200)
nse = 1 - ((mod-ref)**2).sum() / ((ref-ref.mean())**2).sum()
assert abs(1 - (1 - ((ref-ref)**2).sum() / ((ref-ref.mean())**2).sum())) < 1e-15
print(f"3. NSE self=1 OK; noisy model NSE={nse:.3f} (expect ~0.75)")

# ---- 4. effective correlation: AR(1) Neff + betainc p vs scipy
T = 240; phi = 0.7
def ar1(T, phi):
    x = np.zeros(T)
    for t in range(1, T): x[t] = phi * x[t-1] + rng.normal()
    return x
P = []
for _ in range(400):
    x = ar1(T, phi); y = ar1(T, phi)     # independent -> p uniform w/ Neff
    r = np.corrcoef(x, y)[0, 1]
    r1 = np.corrcoef(x[:-1], x[1:])[0, 1]; r2 = np.corrcoef(y[:-1], y[1:])[0, 1]
    neff = np.clip(T * (1 - r1*r2) / (1 + r1*r2), 4, T)
    t = r * np.sqrt((neff - 2) / (1 - r**2))
    from scipy.special import betainc as bi
    p = bi((neff - 2) / 2, 0.5, (neff - 2) / ((neff - 2) + t**2))
    P.append(p)
P = np.array(P)
frac = (P < 0.05).mean()
# naive (uncorrected) for contrast
Pn = []
for _ in range(400):
    x = ar1(T, phi); y = ar1(T, phi)
    r, p = stats.pearsonr(x, y); Pn.append(p)
fn = (np.array(Pn) < 0.05).mean()
print(f"4. effectiveCorr false-positive @5%: corrected {frac:.3f} vs naive {fn:.3f} (expect ~0.05 vs >>0.05)")
assert frac < 0.12 and fn > 0.2

# ---- 5. three-cornered hat (pairwise LS)
T = 20000; sig = np.array([1.0, 2.0, 3.0, 1.5])
common = np.cumsum(rng.normal(size=T)) * 0.05
X = np.array([common + s * rng.normal(size=T) for s in sig]).T   # T x N
N = X.shape[1]
pairs = [(i, j) for i in range(N) for j in range(i+1, N)]
V = np.array([np.var(X[:, i] - X[:, j], ddof=1) for i, j in pairs])
Amat = np.zeros((len(pairs), N))
for k, (i, j) in enumerate(pairs):
    Amat[k, i] = 1; Amat[k, j] = 1
est = np.linalg.lstsq(Amat, V, rcond=None)[0]
est = np.sqrt(np.clip(est, 0, None))
print("5. TCH recovery:", np.round(est, 3), "true", sig)
assert np.allclose(est, sig, rtol=0.08)

# ---- 6. annual amplitude/phase via LS
t = 2019 + np.arange(48) / 12
ytrue = 3.0 * np.cos(2*np.pi*(t - 2019.2)) + 0.1 * rng.normal(size=48)
Am = np.c_[np.ones_like(t), t - t.mean(), np.cos(2*np.pi*t), np.sin(2*np.pi*t)]
coef = np.linalg.lstsq(Am, ytrue, rcond=None)[0]
amp = np.hypot(coef[2], coef[3]); ph = np.arctan2(coef[3], coef[2]) / (2*np.pi)
# y = A cos(2pi(t - t0)) => cos coef = A cos(2pi t0), sin coef = A sin(2pi t0)
print(f"6. annual fit: amp={amp:.3f} (3.0), phase t0={ph:.4f} yr ({0.2:.4f}) "
      f"-> lag days ok: {abs((ph - 0.2)) * 365.25 < 2}")
assert abs(amp - 3.0) < 0.1 and abs(ph - 0.2) * 365.25 < 2
print("ALL COMPARISON NUMERICS VALIDATED")
