import numpy as np

# Cross-check the parity split + polynomial removal logic (same algorithm as
# shDestripe.m local functions) independently in Python.
nmax = 40
m = 10
n = np.arange(m, nmax+1)
signal_smooth = 0.01*(n-20.0)**2 * 1e-9

even = n[n % 2 == 0]
odd  = n[n % 2 == 1]

def poly_remove(nsub, y, order):
    coeff = np.polyfit(nsub, y, order)
    trend = np.polyval(coeff, nsub)
    return y - trend

y_even = signal_smooth[n % 2 == 0]
y_odd  = signal_smooth[n % 2 == 1]
res_even = poly_remove(even, y_even, 2)
res_odd  = poly_remove(odd, y_odd, 2)

print("max |residual| even:", np.max(np.abs(res_even)))
print("max |residual| odd :", np.max(np.abs(res_odd)))
assert np.max(np.abs(res_even)) < 1e-13
assert np.max(np.abs(res_odd)) < 1e-13
print("Destriping logic (parity split + polynomial removal) validated: quadratic trend fully absorbed by order-2 fit.")
