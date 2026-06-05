#!/usr/bin/env python3
"""
Spectral gap of the desingularized stability operator (item 1, the crux of
Conjecture 2.4): is the profile ASYMPTOTICALLY stable (essential spectrum strictly
left of 0) or only MARGINALLY stable?

The transport -2 d_z is CONVECTIVE (carries perturbations to +z, the far field, at
speed 2).  In the plain norm L^2(dY) (= G in L^2(dz)) the essential spectrum is on
the imaginary axis -> MARGINAL, no absolute decay (perturbations propagate out
without decaying pointwise; the discrete operator shows Re -> 0^-).  A finite-matrix
weight cannot move this (similarity preserves eigenvalues).

The gap lives in a WEIGHTED (convective) norm.  Analytic essential spectrum in
L^2(e^{-2 gamma z} dz) (edge decay required; far-field outflow allowed -- the only
admissible direction for an outflow operator): the symbol of -2 d_z conjugated by
e^{-gamma z} is -2 d_z - 2 gamma, so

        essential spectrum:  Re = -2 gamma  <  0   for gamma > 0   =>  GAP = 2 gamma.

Physical/numerical confirmation here = the SEMIGROUP DECAY RATE, which equals the
weighted spectral abscissa.  Prediction: a localized perturbation's UNWEIGHTED norm
stays ~flat (marginal), while its EDGE-WEIGHTED norm ||e^{-gamma z} G|| decays at
rate exactly 2 gamma.  (No unstable point mode exists -- established earlier -- so
the abscissa is set by the essential spectrum.)
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from stability_edge_z import build_blocks_z  # noqa: E402


def evolve(M, z, G0, taus):
    """G(tau) = expm(M tau) G0 via eigendecomposition (linear, exact)."""
    vals, V = np.linalg.eig(M)
    c = np.linalg.solve(V, G0)
    out = []
    for t in taus:
        out.append(V @ (np.exp(vals * t) * c))
    return out


def norms(z, dz, G, gammas):
    w = {g: float(np.sqrt(np.sum(np.exp(-2 * g * z) * np.abs(G) ** 2) * dz))
         for g in gammas}
    w[0.0] = float(np.sqrt(np.sum(np.abs(G) ** 2) * dz))   # unweighted L2(dz)
    return w


def main():
    print("=" * 78)
    print("  SPECTRAL GAP via convective semigroup decay (edge-resolved operator)")
    print("=" * 78)
    B = build_blocks_z(N=1400, zmin=-14.0, zmax=12.0, edge_bc=True)
    M, z = B["M_GG"], B["z"]
    dz = z[1] - z[0]

    # localized bulk shape perturbation (well away from both ends)
    G0 = np.exp(-((z - 1.0) / 1.0) ** 2)

    gammas = [0.25, 0.5, 0.75]
    taus = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    Gt = evolve(M, z, G0, taus)
    series = {g: [] for g in [0.0] + gammas}
    for G in Gt:
        nm = norms(z, dz, G, gammas)
        for g in series:
            series[g].append(nm[g])

    print("\n  norm ratios ||.||(tau)/||.||(0):")
    hdr = "   tau   " + "  unweighted" + "".join(f"   gamma={g:.2f}" for g in gammas)
    print(hdr)
    for i, t in enumerate(taus):
        row = f"   {t:.1f}   {series[0.0][i]/series[0.0][0]:11.4f}"
        for g in gammas:
            row += f"   {series[g][i]/series[g][0]:10.4f}"
        print(row)

    # fit decay rate of each weighted norm over the clean window (before bump exits)
    print("\n  fitted decay rate  -d/dtau ln||.||  vs predicted gap 2*gamma:")
    win = (taus >= 0.5) & (taus <= 2.5)
    tt = taus[win]
    for g in [0.0] + gammas:
        y = np.log(np.array(series[g])[win])
        rate = -np.polyfit(tt, y, 1)[0]
        pred = 2 * g
        tag = "(marginal)" if g == 0.0 else f"(predict 2*gamma={pred:.2f})"
        print(f"   gamma={g:.2f}:  measured decay rate = {rate:+.3f}   {tag}")

    print("\n" + "=" * 78)
    print("  READ: unweighted norm ~ flat (marginal, no absolute decay); edge-weighted")
    print("  norms decay at rate ~ 2*gamma -> the operator has a SPECTRAL GAP of 2*gamma")
    print("  in L^2(e^{-2 gamma z} dz).  Asymptotic (convective, weighted) stability of")
    print("  the Omega-shape sector, modulo the symmetry modes.  Consistent with")
    print("  Conjecture 2.4; this is the standard weighted-norm sense for transport-")
    print("  dominated self-similar operators (NOT absolute decay in L^2(dY)).")
    print("=" * 78)


if __name__ == "__main__":
    main()
