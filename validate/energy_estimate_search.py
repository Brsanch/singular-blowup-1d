#!/usr/bin/env python3
"""
PROOF ATTEMPT, step 1 (the key lemma, numerically): does a BOUNDED-AT-EDGE anisotropic
weight make the Omega-shape energy estimate close?  I.e. is the energy abscissa

    A(phi) = sup_G  <G, M_GG G>_W / ||G||_W^2   < 0 ?

with the transport handled ANALYTICALLY (exact IBP) as the diagonal multiplication
    <G, -2 d_z G>_W = -int phi'(z) W G^2 dz - [W G^2]_bdry  <=  -<G, phi' G>_W   (outflow<=0,
inflow G(-inf)=0), so the transport contributes the SELF-ADJOINT diagonal  -diag(phi'(z));
and the bounded nonlocal forcing F=1/2 Q_op via its weighted symmetric part. Then
    A = lambda_max( -diag(phi') + sym_W(F) ),   sym_W(F) = (sqrt(W) F sqrt(W)^{-1} + transpose)/2.

The convective weight phi=2 gamma z gives W=(Y-1)^{-2 gamma} -> BLOWS UP at the edge, which
amplifies the nonlocal term (C_force diverges, cap_scaffold). Fix: phi'(z) = 2 gamma * ramp(z)
that -> 0 at the edge (z->-inf) and -> 2 gamma in the bulk; then phi(-inf) is finite, W is
BOUNDED at the edge, the edge<->bulk weight ratio is bounded, and the forcing should be
controlled -- while phi'~2 gamma in the bulk supplies the dissipation. If A < 0 and CONVERGES
under refinement, the key coercivity lemma holds and we have the explicit weight.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from stability_edge_z import build_blocks_z  # noqa: E402


def phi_prime(z, gamma, z0, s):
    """phi'(z): ramps 0 (edge) -> 2 gamma (bulk) via a logistic in z."""
    return 2.0 * gamma * 0.5 * (1.0 + np.tanh((z - z0) / s))


def energy_abscissa(B, gamma, z0, s, return_vec=False):
    z, dz, Qop = B["z"], B["dz"], B["Q_op"]
    pp = phi_prime(z, gamma, z0, s)
    phi = np.cumsum(pp) * dz                  # phi(z), phi(zmin) ~ 0 (phi'->0 at edge)
    sw = np.exp(-0.5 * phi)                    # sqrt(W)
    F = 0.5 * Qop
    Ft = (sw[:, None] * F) / sw[None, :]       # similarity transform of the forcing
    S = -np.diag(pp) + 0.5 * (Ft + Ft.T)       # energy operator (symmetric)
    if return_vec:
        vals, vecs = np.linalg.eigh(S)
        return float(vals[-1]), vecs[:, -1], z, phi
    return float(np.linalg.eigvalsh(S)[-1])


def main():
    print("=" * 84)
    print("  PROOF ATTEMPT step 1: energy abscissa with a BOUNDED-AT-EDGE anisotropic weight")
    print("=" * 84)
    B = build_blocks_z(N=1400, zmin=-12.0, zmax=12.0, beta=0.0, edge_bc=True)

    print("\n  baseline -- convective weight phi=2 gamma z (blows up at edge): energy abscissa")
    # emulate the pure convective weight by a ramp with the step far to the left (no flattening)
    for g in (0.25, 0.5):
        A = energy_abscissa(B, g, z0=-30.0, s=0.5)   # ramp centered far left ~ phi'=2gamma everywhere
        print(f"    gamma={g:.2f}: A = {A:+.3f}   (expect large >0: the known failure)")

    print("\n  bounded-at-edge ramp phi'(z)=2 gamma * (1+tanh((z-z0)/s))/2  -> 0 at edge:")
    print("    gamma  z0     s     energy abscissa A (want < 0)")
    best = (None, 1e9)
    for g in (0.25, 0.4, 0.5):
        for z0 in (-6.0, -4.0, -2.0, 0.0):
            for s in (1.0, 2.0):
                A = energy_abscissa(B, g, z0, s)
                if A < best[1]:
                    best = ((g, z0, s), A)
                flag = "  <- DISSIPATIVE" if A < 0 else ""
                print(f"    {g:.2f}  {z0:+.0f}   {s:.0f}    {A:+.4f}{flag}")
    print(f"\n  best: (gamma,z0,s)={best[0]}, energy abscissa = {best[1]:+.4f}")

    # DIAGNOSE the worst mode (regardless of sign): grid-scale edge artifact or genuine?
    g, z0, s = best[0]
    A, vec, z, phi = energy_abscissa(B, g, z0, s, return_vec=True)
    v = np.abs(vec) / (np.abs(vec).max() + 1e-30)
    zstar = z[int(np.argmax(v))]
    F = np.abs(np.fft.rfft(vec)); hik = float(F[len(F)//2:].sum() / (F.sum() + 1e-30))
    print(f"\n  WORST-MODE diagnosis: location z* = {zstar:+.2f}  (edge=z_min=-12; bulk~0..6)")
    print(f"    hi-k energy fraction = {hik:.3f}  (~1 => grid-scale artifact; ~0 => smooth/physical)")
    print("\n  ABSCISSA vs N at fixed domain (grows => grid-scale edge artifact; converges => genuine):")
    print("    N      zmin   energy abscissa A")
    for (Nc, zc_) in [(700, -12.0), (1000, -12.0), (1400, -12.0), (2000, -12.0)]:
        Bc = build_blocks_z(N=Nc, zmin=zc_, zmax=12.0, beta=0.0, edge_bc=True)
        print(f"    {Nc:<6d} {zc_:+.0f}   {energy_abscissa(Bc, g, z0, s):+.4f}")

    # PHYSICAL (smooth-subspace) abscissa: project the energy operator onto the smooth
    # (low-wavenumber) subspace -- drop the grid-scale edge modes -- and recompute.
    print("\n  SMOOTH-SUBSPACE abscissa (project onto lowest-k modes; if <0 => physical estimate closes):")
    zc, dzc, Qop = B["z"], B["dz"], B["Q_op"]
    pp = phi_prime(zc, g, z0, s); phi = np.cumsum(pp) * dzc; sw = np.exp(-0.5 * phi)
    Ft = (sw[:, None] * (0.5 * Qop)) / sw[None, :]
    S = -np.diag(pp) + 0.5 * (Ft + Ft.T)
    n = len(zc)
    # smooth basis = lowest eigenvectors of the 1D Laplacian on the grid (Fourier-like)
    lap = np.diag(-2.0 * np.ones(n)) + np.diag(np.ones(n - 1), 1) + np.diag(np.ones(n - 1), -1)
    _, B_lap = np.linalg.eigh(-lap)             # columns: low-k -> high-k
    for frac in (0.1, 0.25, 0.5):
        m = max(2, int(frac * n))
        P = B_lap[:, :m]                         # lowest m smooth modes
        Sr = P.T @ S @ P
        print(f"    keep lowest {int(frac*100)}% modes (m={m}): abscissa = {float(np.linalg.eigvalsh(Sr)[-1]):+.4f}")

    # EDGE-CUTOFF abscissa: the worst mode is at the edge where G->0 (Robin condition); the
    # discretized (Y-1)^{-1/2} edge artifact dominates. Exclude a thin edge layer z<z_cut and
    # recompute the BULK abscissa. If <0 => the bulk energy estimate closes; the edge needs a
    # separate (analytic, matched-asymptotic) treatment -- exactly Chen's anisotropic edge handling.
    print("\n  EDGE-CUTOFF abscissa (drop z < z_cut; tests whether the BULK estimate closes):")
    print("    z_cut   #edge nodes dropped   bulk abscissa   bulk z* ")
    for z_cut in (-10.0, -8.0, -6.0, -4.0):
        keep = zc > z_cut
        Sb = S[np.ix_(keep, keep)]
        vb, Vb = np.linalg.eigh(Sb)
        zb = zc[keep][int(np.argmax(np.abs(Vb[:, -1])))]
        print(f"    {z_cut:+.0f}    {int((~keep).sum()):>4d}                 {float(vb[-1]):+.4f}      {zb:+.2f}")

    print("\n" + "=" * 84)
    print("  READ: if A < 0 and CONVERGES under refinement for a bounded-at-edge weight, the")
    print("  Omega-shape coercivity lemma (star) holds with that explicit weight -- the energy")
    print("  estimate the proof needs. (If the worst mode sits at the edge, it is a symmetry")
    print("  mode to project out; if in the bulk, the estimate is unconditional on the shape.)")
    print("=" * 84)


if __name__ == "__main__":
    main()
