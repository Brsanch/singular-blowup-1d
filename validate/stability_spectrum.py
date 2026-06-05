#!/usr/bin/env python3
"""
Convergence study for the leading eigenvalue of the singular-profile stability
operator.  The question: is the unstable lambda ~ +0.14 (frozen, mu=nu=0) a
GENUINE instability or a regularization/mesh artifact?

For each (N, delta, box) we report the top real eigenvalues, the leading
unstable lambda (excluding the ~0 symmetry mode), and where its eigenfunction
concentrates (edge X=1 vs bulk) + whether it is grid-resolved.
"""
import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from stability import build_operator, spectrum  # noqa: E402


def leading(L, X):
    vals, vecs = spectrum(L)
    N = len(X)
    # exclude near-zero symmetry mode(s) when reporting the "instability"
    reals = vals.real
    # leading eigenvalue overall and leading EXCLUDING |lambda|<3e-3 (symmetry)
    lam_top = vals[0]
    nontriv = [i for i in range(len(vals)) if abs(vals[i]) > 3e-3]
    i_unstable = nontriv[0] if nontriv else 0
    lam_u = vals[i_unstable]
    v = vecs[:N, i_unstable]                 # omega-component of eigenvector
    vmag = np.abs(v)
    xloc = X[int(np.argmax(vmag))]
    # resolution of the mode: fraction of |v| in the top 1/6 wavenumbers
    F = np.abs(np.fft.fft(v))
    k = np.fft.fftfreq(N, d=1.0 / N)
    hi = np.abs(k) > N // 6
    tail = float(np.sum(F[hi] ** 2) / (np.sum(F ** 2) + 1e-30))
    return lam_top, lam_u, xloc, tail, vals[:5]


def main():
    print("=" * 78)
    print("  STABILITY SPECTRUM CONVERGENCE  (frozen mu=nu=0)")
    print("=" * 78)

    print("\n  (A) delta-dependence at fixed N, box [-8,40]:")
    print("   delta     lambda_max(nontriv)   mode x*    hi-k frac")
    for delta in (0.14, 0.10, 0.07, 0.05, 0.035):
        L, X = build_operator(N=1100, x0=-8.0, x1=40.0, delta=delta)
        ltop, lu, xloc, tail, _ = leading(L, X)
        print(f"   {delta:.3f}    {lu.real:+.5f} {lu.imag:+.5f}i    "
              f"{xloc:+.3f}    {tail:.2e}", flush=True)

    print("\n  (B) mesh-dependence at fixed delta=0.06, box [-8,40]:")
    print("   N        lambda_max(nontriv)   mode x*    hi-k frac")
    for N in (700, 1000, 1400, 1900):
        L, X = build_operator(N=N, x0=-8.0, x1=40.0, delta=0.06)
        ltop, lu, xloc, tail, _ = leading(L, X)
        print(f"   {N:5d}    {lu.real:+.5f} {lu.imag:+.5f}i    "
              f"{xloc:+.3f}    {tail:.2e}", flush=True)

    print("\n  (C) box-dependence at fixed N-density, delta=0.06:")
    print("   box            lambda_max(nontriv)   mode x*")
    for (x0, x1) in ((-6.0, 25.0), (-8.0, 40.0), (-12.0, 70.0)):
        dens = 1100 / 48.0
        N = int(dens * (x1 - x0)) | 1
        L, X = build_operator(N=N, x0=x0, x1=x1, delta=0.06)
        ltop, lu, xloc, tail, _ = leading(L, X)
        print(f"   [{x0:.0f},{x1:.0f}] N={N:4d}   {lu.real:+.5f} {lu.imag:+.5f}i"
              f"    {xloc:+.3f}", flush=True)

    print("\n" + "=" * 78)
    print("  READ: if lambda_max(nontriv) drifts with delta -> regularization")
    print("  artifact (free-boundary mode). If it converges in N AND delta and")
    print("  the eigenmode is bulk + low hi-k frac -> a genuine instability.")
    print("=" * 78)


if __name__ == "__main__":
    main()
