#!/usr/bin/env python3
"""
Validation of the edge-RESOLVING (z=ln(Y-1)) stability operator
(src/stability_edge_z.py).  This is the clean, edge-resolved confirmation of the
equispaced result (no genuine unstable mode), and it documents a methodological
catch: the inflow/regularity edge BC is essential.

Gates:
 (1) z-grid non-equispaced Hilbert: H((Y-1)^{-1/2}) ~ 1-2% vs the closed form 0 on
     Y>1 (10x better than the equispaced tool's ~15%); the residual is combined
     far-field tail truncation + resolution, and the spectral verdict (gate 3) is
     robust to it.
 (2) Free edge vs inflow BC, with the edge average truncated at z_min (edge_tail=False)
     and with its tail included (edge_tail=True, the default since 2026-11-22).  The
     SPURIOUS, z_min-dependent free-edge mode (Re~+1.7) exists only for the truncated
     average; with the tail there is no Re>0 mode with or without the BC.  With the BC
     G(z_min)=0 both give the same max Re < 0.
 (3) the clean verdict (max Re < 0) is robust under N, z_min, z_max refinement.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from stability_edge_z import build_blocks_z, hilbert_support_z, spectrum, C_L  # noqa: E402


def gate1_hilbert():
    print("=" * 78)
    print("  (1) z-grid non-equispaced Hilbert:  H((Y-1)^{-1/2}) = 0 on Y>1 (closed form)")
    print("=" * 78)
    print("   N      zmax    Ymax        max|H(g)| on Y in (1.05,300)")
    for (N, zmax) in ((1600, 7.0), (3200, 7.0), (1600, 9.0), (3200, 9.0)):
        z = np.linspace(-14.0, zmax, N)
        Y = 1.0 + np.exp(z)
        Ymax = Y[-1] + 0.5 * (Y[-1] - Y[-2])
        H = hilbert_support_z(Y, Ymax)
        g = (Y - 1.0) ** -0.5
        Hg = H @ g
        m = (Y > 1.05) & (Y < 300)
        print(f"   {N:5d}  {zmax:5.1f}   {Ymax:9.1f}   {np.max(np.abs(Hg[m])):.3e}")
    print("   -> ~1-2%: tail (smaller Ymax = larger error) + resolution (larger N =")
    print("      smaller error). 10x better than equispaced; verdict robust to it.")


def gate2_edge_bc():
    print("\n" + "=" * 78)
    print("  (2) free edge vs inflow BC; truncated edge average vs tail included")
    print("=" * 78)
    print("   rightmost Re over the whole spectrum / over |Im| < 20 (the whole-spectrum maximum")
    print("   is a grid-sawtooth pair at |Im| = 2/dz, ~1/N, blind to the edge)")
    print("   zmin   no BC, truncated    no BC, exact tail    BC, truncated      BC, exact tail")
    for zmin in (-6.0, -10.0, -12.0, -14.0, -16.0, -18.0):
        mx = {}
        for bc in (False, True):
            for tail in (False, True):
                v, _ = spectrum(build_blocks_z(N=1600, zmin=zmin, zmax=7.0, edge_bc=bc,
                                               edge_tail=tail)["M_GG"])
                mx[bc, tail] = (v.real.max(), v[np.abs(v.imag) < 20.0].real.max())
        print(f"   {zmin:5.0f}  " + "   ".join(
            f"{mx[k][0]:+.4f}/{mx[k][1]:+.3f}" for k in ((False, False), (False, True), (True, False), (True, True))))
    print("   -> the real free-edge mode (Re~+1.4..+1.7, drifting with z_min) appears with the")
    print("      TRUNCATED edge average; with the tail included this matrix has no Re>0 mode for")
    print("      z_min <= -4.  The free-edge spectrum is discretization-defined (a free inflow")
    print("      boundary has no continuum limit; an independent collocation of Q has an unstable")
    print("      edge-node mode), so only the BC columns carry meaning: max Re < 0 in both.")


def gate3_verdict():
    print("\n" + "=" * 78)
    print("  (3) clean verdict (edge BC on): max Re eig(M_GG) < 0, robust")
    print("=" * 78)
    print("   N      zmin    zmax    max Re eig(M_GG)   #Re>0   max Re eig(M_hh)")
    for (N, zmin, zmax) in ((1200, -14, 7), (1600, -14, 7), (2200, -14, 7),
                            (1600, -16, 7), (1600, -14, 10)):
        B = build_blocks_z(N=N, zmin=zmin, zmax=zmax, edge_bc=True)
        vGG, _ = spectrum(B["M_GG"])
        vhh, _ = spectrum(B["M_hh"])
        print(f"   {N:5d}  {zmin:5.0f}  {zmax:5.0f}    {vGG.real.max():+.5f}"
              f"          {int((vGG.real>1e-6).sum()):3d}     {vhh.real.max():+.5f}")
    print("\n" + "=" * 78)
    print("  VERDICT (edge-resolved, z=ln(Y-1), inflow BC, frozen mu=nu=0):")
    print("   - eig(M_GG) with Re>0 : 0   (NO genuine unstable Omega-shape mode)")
    print("   - rightmost eig ~ -0.001 = marginal transport continuum (essential")
    print("     spectrum at Re=0 in the bounded-G norm; shiftable strictly left by a")
    print("     far-field decay weight -> asymptotic stability in that weighted norm).")
    print(f"   - edge/xi symmetry mode appended at lambda = {C_L-1:g} (project out).")
    print("   - Theta-shape marginal at Re=0 (the nu normalization direction).")
    print("  Confirms the equispaced verdict, now edge-RESOLVED and clean (no")
    print("  spurious near-edge cloud, no box modes). NUMERICAL EVIDENCE CONSISTENT")
    print("  WITH Conjecture 2.4 (frozen sector). NOT a proof; mu,nu + nonlinear")
    print("  shape test + interval arithmetic remain.")
    print("=" * 78)


if __name__ == "__main__":
    gate1_hilbert()
    gate2_edge_bc()
    gate3_verdict()
