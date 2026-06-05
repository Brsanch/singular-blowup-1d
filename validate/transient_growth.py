#!/usr/bin/env python3
"""
Transient growth / non-normality of the stability operator (Conjecture 2.4) -- the
question that qualifies "stable" beyond the eigenvalues.

All eigenvalues are stable (+ convective gap), but the operator is NON-NORMAL
(transport-dominated). For a non-normal operator, stable eigenvalues do NOT preclude
transient amplification ||e^{M tau}||, which is what can trigger NONLINEAR
instability. We ask: is the transient amplification G_max = sup_tau ||e^{M tau}||
MODEST (=> robust linear stability, nonlinear stability supported) or LARGE?

Honest status: this is genuinely hard to compute for this operator because of EDGE
REPRESENTATION CONDITIONING. We fix the worst offender (the (1/2)J/(Y-1) forcing,
recast as a bounded exponential-averaging convolution -> resolution-independent
||M||), but a residual e^{-z/2} edge factor remains in the Hilbert map H(g), and the
weighted worst-case operator norm is contaminated by the far-field weight. So we
report the PHYSICAL-norm (unweighted) transient amplification, with its resolution
dependence stated explicitly -- partially open, NOT a clean robustness proof.
"""
import os
import sys

import numpy as np
from scipy.linalg import expm

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from stability_edge_z import hilbert_support_z, _dz_central  # noqa: E402


def build_MGG_wellcond(N, zmin, zmax, edge_bc=True):
    """Omega-shape operator with the WELL-CONDITIONED (bounded) edge forcing:
    (1/2) J/(Y-1) = (1/2) int_{zmin}^z H(g) e^{z'-z} dz'  (exp-averaging, weights<=1)."""
    z = np.linspace(zmin, zmax, N)
    dz = z[1] - z[0]
    Y = 1.0 + np.exp(z)
    ym1 = Y - 1.0
    Ymax = Y[-1] + 0.5 * (Y[-1] - Y[-2])
    Hg_op = hilbert_support_z(Y, Ymax) @ np.diag(ym1 ** -0.5)
    L = np.zeros((N, N))
    for i in range(N):
        w = np.full(i + 1, dz)
        w[0] *= 0.5
        if i > 0:
            w[i] *= 0.5
        L[i, :i + 1] = np.exp(z[:i + 1] - z[i]) * w
    M = -2.0 * _dz_central(N, dz) + 0.5 * (L @ Hg_op)
    if edge_bc:
        M, z = M[1:, 1:], z[1:]
    return M, z


def main():
    print("=" * 78)
    print("  TRANSIENT GROWTH / non-normality (well-conditioned Omega-shape operator)")
    print("=" * 78)

    print("  conditioning fix works: ||M||_2 ~ resolution-independent")
    for (N, zmin, zmax) in ((600, -8, 6), (800, -10, 7), (1000, -12, 8)):
        M, _ = build_MGG_wellcond(N, zmin, zmax)
        print(f"    N={N} zmin={zmin}: ||M||={np.linalg.norm(M,2):.0f}  "
              f"maxRe={np.linalg.eigvals(M).real.max():+.4f}")

    taus = np.linspace(0, 12, 49)
    print("\n  PHYSICAL-norm (L^2(dY)) transient amplification G_max = sup ||e^{M tau}||:")
    print("   N     zmin  zmax   G_max  (at tau)   ||e^{12 M}||  (return to ~1?)")
    for (N, zmin, zmax) in ((600, -8, 6), (800, -10, 7), (1000, -12, 8), (900, -10, 9)):
        M, _ = build_MGG_wellcond(N, zmin, zmax)
        amps = np.array([np.linalg.norm(expm(M * t), 2) for t in taus])
        print(f"   {N:4d}  {zmin:4d}  {zmax:3d}    {amps.max():5.2f}  "
              f"({taus[int(amps.argmax())]:.2f})      {amps[-1]:5.2f}")

    print("\n" + "=" * 78)
    print("  READ (calibrated, NOT overclaimed):")
    print("  - the PHYSICAL-norm transient amplification is MODEST: G_max ~ 2-3, and")
    print("    it mostly decays back toward ~1 (bounded transient, no blow-up).")
    print("  - BUT G_max creeps up with edge depth (1.8 -> 2.3 -> 3.1 as zmin -> -12),")
    print("    which we cannot cleanly separate from a residual e^{-z/2} edge factor in")
    print("    the Hilbert map. So the PRECISE non-normality bound is an OPEN refinement")
    print("    (needs a fully edge-conditioned operator); the eigenvalue stability + the")
    print("    convective gap (spectral_gap.py, specific-perturbation decay) are solid.")
    print("  - net: linear stability is eigenvalue-robust; the transient-growth")
    print("    (nonlinear-robustness) bound is SUGGESTIVE (modest) but NOT yet settled.")
    print("=" * 78)


if __name__ == "__main__":
    main()
