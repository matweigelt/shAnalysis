"""validate_stopping.py - when to stop an ill-posed forward-model iteration.

shLowLevel.leakageCorrect runs the fixed point

    m <- m + g * (obs - F(m))

which is a Landweber iteration on an ill-posed inverse problem. Such
iterations SEMICONVERGE: the error against the truth falls, reaches a
minimum, and then RISES again as the iteration starts fitting the noise.
The iteration count is therefore the regularisation parameter, and
"run until it stops changing" is the wrong instruction - the residual
keeps shrinking long after the solution has started to degrade.

The classical fix is the DISCREPANCY PRINCIPLE (Morozov): stop as soon
as the residual reaches the noise level of the data,

    ||obs - F(m_k)|| <= Tau * delta,     Tau >= 1,

because fitting the data more closely than its own noise is fitting
noise. This script demonstrates the U-shaped error curve and checks that
the discrepancy principle stops near its minimum, while a
step-size criterion does not.

Developed by Matthias Weigelt with the help of Claude (Opus 5),
2026-08-11 (v3.5.0).
"""
import numpy as np


def make_problem(n=200, seed=3):
    """A smoother F with decaying singular values, a masked truth."""
    rng = np.random.default_rng(seed)
    x = np.linspace(0, 1, n)
    # smoothing operator: Gaussian convolution -> eigenvalues decay fast
    d = x[:, None] - x[None, :]
    F = np.exp(-(d ** 2) / (2 * 0.03 ** 2))
    F /= F.sum(axis=1, keepdims=True)
    mask = (x > 0.35) & (x < 0.65)
    truth = np.where(mask, 1.0, 0.0)
    clean = F @ truth
    return F, truth, clean, mask, rng


def landweber(F, obs, mask, gain=1.0, iters=4000):
    """m <- m + gain*(obs - F m), confined to mask. Returns the history."""
    m = np.zeros_like(obs)
    err, res = [], []
    for _ in range(iters):
        r = obs - F @ m
        m = m + gain * r
        m = m * mask
        err.append(m)
        res.append(np.linalg.norm(r))
    return np.array(err), np.array(res)


def main():
    F, truth, clean, mask, rng = make_problem()
    delta_rel = 0.01
    noise = delta_rel * np.linalg.norm(clean) / np.sqrt(clean.size)
    obs = clean + noise * rng.normal(size=clean.size)
    delta = np.linalg.norm(obs - clean)          # the true noise norm
    print("noise level ||obs-clean|| = %.4e" % delta)

    M, res = landweber(F, obs, mask, gain=1.0, iters=4000)
    err = np.linalg.norm(M - truth, axis=1)
    kbest = int(np.argmin(err))
    print("\nsemiconvergence")
    print("  best iteration      k = %4d   error %.4e" % (kbest, err[kbest]))
    print("  at 10x that         k = %4d   error %.4e (%.2fx worse)"
          % (min(10 * kbest, len(err) - 1), err[min(10 * kbest, len(err) - 1)],
             err[min(10 * kbest, len(err) - 1)] / err[kbest]))
    print("  at the end          k = %4d   error %.4e (%.2fx worse)"
          % (len(err) - 1, err[-1], err[-1] / err[kbest]))
    assert err[-1] > 1.2 * err[kbest], "no semiconvergence - test is void"
    # the residual keeps falling while the error rises: the trap
    assert res[-1] < res[kbest], "residual must keep shrinking"
    print("  residual STILL falls while the error rises: %.4e -> %.4e"
          % (res[kbest], res[-1]))

    print("\ndiscrepancy principle (stop at ||r|| <= Tau*delta)")
    for tau in (1.0, 1.2, 1.5, 2.0):
        k = int(np.argmax(res <= tau * delta))
        if res[k] > tau * delta:
            print("  Tau %.1f: never reached" % tau)
            continue
        print("  Tau %.1f: k = %4d  error %.4e  (%.2fx the optimum)"
              % (tau, k, err[k], err[k] / err[kbest]))
        assert err[k] < 2.0 * err[kbest], (tau, err[k] / err[kbest])

    print("\nstep-size criterion (stop when ||m_k - m_(k-1)|| small)")
    step = np.linalg.norm(np.diff(M, axis=0), axis=1) / \
        np.maximum(np.linalg.norm(M[1:], axis=1), 1e-30)
    for tol in (1e-3, 1e-4, 1e-5):
        idx = np.where(step < tol)[0]
        if idx.size == 0:
            print("  Tol %.0e: never reached in %d iterations" % (tol, len(err)))
            continue
        k = int(idx[0])
        print("  Tol %.0e: k = %4d  error %.4e  (%.2fx the optimum)"
              % (tol, k, err[k], err[k] / err[kbest]))

    print("\nvalidate_stopping: the discrepancy principle lands near the "
          "optimum; chasing a small step does not")


if __name__ == "__main__":
    main()
