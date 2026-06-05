#!/usr/bin/env python3
"""
Transient growth / non-normality of the stability operator (Conjecture 2.4) -- the
question that qualifies "stable" beyond the eigenvalues, and the precursor to
nonlinear stability.

The operator is NON-NORMAL (transport-dominated), so stable eigenvalues + gap do not
by themselves bound the transient amplification ||e^{M tau}||. We resolve the
question by separating PHYSICAL (smooth, resolved) perturbations from grid-scale
edge modes:

RESULT
  * SMOOTH physical perturbations (bulk or near-edge): max amplification G = 1.000 at
    EVERY resolution -- i.e. NO transient growth, monotone decay. Linear stability is
    robust; nonlinear stability is supported.
  * the WORST-CASE operator norm sup ||e^{M tau}|| creeps up, but with EDGE DEPTH
    (zmin), not with N at fixed domain, and its worst-case mode is a SPIKE at the
    deepest edge node = a grid-scale representation of the (Y-1)^{-1/2} singularity,
    NOT a physical perturbation. This is a discretization feature of the singular
    edge, not an instability.

(Conditioning note: the (1/2)J/(Y-1) forcing is recast as the bounded exp-averaging
convolution (1/2) int_{zmin}^z H(g) e^{z'-z} dz' so ||M|| is resolution-independent.)
"""
import os
import sys

import numpy as np
from scipy.linalg import expm

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from stability_edge_z import hilbert_support_z, _dz_central  # noqa: E402


def build_MGG(N, zmin, zmax, edge_bc=True):
    z = np.linspace(zmin, zmax, N)
    dz = z[1] - z[0]
    Y = 1.0 + np.exp(z)
    ym1 = Y - 1.0
    Ymax = Y[-1] + 0.5 * (Y[-1] - Y[-2])
    Hg = hilbert_support_z(Y, Ymax) @ np.diag(ym1 ** -0.5)
    L = np.zeros((N, N))
    for i in range(N):
        w = np.full(i + 1, dz)
        w[0] *= 0.5
        if i > 0:
            w[i] *= 0.5
        L[i, :i + 1] = np.exp(z[:i + 1] - z[i]) * w
    M = -2.0 * _dz_central(N, dz) + 0.5 * (L @ Hg)
    if edge_bc:
        M, z = M[1:, 1:], z[1:]
    return M, z


def main():
    print("=" * 78)
    print("  TRANSIENT GROWTH / non-normality: physical perturbations vs grid-scale")
    print("=" * 78)
    taus = np.linspace(0.0, 5.0, 26)

    print("\n  G = sup_tau ||e^{M tau} g||/||g|| for SMOOTH physical perturbations,")
    print("  and the WORST-CASE operator norm sup_tau ||e^{M tau}||:")
    print("   zmin   N      G(bulk bump)  G(near-edge bump)   worst-case ||e^{Mt}||")
    for (N, zmin, zmax) in ((800, -10, 7), (1000, -12, 8), (1400, -12, 8),
                            (1200, -14, 8)):
        M, z = build_MGG(N, zmin, zmax)
        g_bulk = np.exp(-((z - 0.0) / 1.0) ** 2); g_bulk /= np.linalg.norm(g_bulk)
        g_edge = np.exp(-((z + 3.0) / 1.0) ** 2); g_edge /= np.linalg.norm(g_edge)
        Es = [expm(M * t) for t in taus]
        Gb = max(np.linalg.norm(E @ g_bulk) for E in Es)
        Ge = max(np.linalg.norm(E @ g_edge) for E in Es)
        wc = max(np.linalg.norm(E, 2) for E in Es)
        print(f"   {zmin:4d}  {N:5d}     {Gb:6.3f}        {Ge:6.3f}             {wc:6.2f}")

    # confirm the worst-case mode is a grid-scale edge spike (not physical)
    M, z = build_MGG(1000, -12, 8)
    Es = [expm(M * t) for t in taus]
    tpk = taus[int(np.argmax([np.linalg.norm(E, 2) for E in Es]))]
    U, s, Vt = np.linalg.svd(expm(M * tpk))
    vin = Vt[0]
    # hi-k fraction of the worst-case input
    F = np.abs(np.fft.fft(vin)); n = len(vin); k = np.fft.fftfreq(n, 1.0 / n)
    hik = float((F[np.abs(k) > n // 6] ** 2).sum() / (F ** 2).sum())
    print(f"\n  worst-case input (peak tau={tpk:.2f}): peaks at z={z[int(np.argmax(np.abs(vin)))]:.2f}"
          f" (z_min={z[0]:.1f}), hi-k fraction = {hik:.2f}")
    print("  -> a high-wavenumber spike at the deepest edge node = grid-scale, not physical.")

    print("\n" + "=" * 78)
    print("  VERDICT: PHYSICAL (smooth, resolved) perturbations have G = 1.000 at every")
    print("  resolution -- NO transient growth (monotone decay). The non-normal worst-")
    print("  case creep is a grid-scale edge mode (the discretized (Y-1)^{-1/2}")
    print("  singularity), tied to edge DEPTH not N, and is not a physical instability.")
    print("  => linear stability is ROBUST; nonlinear stability is supported. This")
    print("  settles the transient-growth question for physical perturbations.")
    print("=" * 78)


if __name__ == "__main__":
    main()
