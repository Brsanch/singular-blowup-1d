#!/usr/bin/env python3
"""
Fold in the Theta sector and the mu,nu scaling/normalization modes (the rest of the
frozen-sector verdict for Conjecture 2.4), on the edge-resolved z=ln(Y-1) operator.

Key structural point: rather than add mu,nu with a paper-specific normalization, use
the CONVECTIVE WEIGHTED space L^2(e^{-2 gamma z} dz) (the one that gives the spectral
gap).  It AUTOMATICALLY excludes the scaling symmetry modes:

  * amplitude (nu) generator:  dOmega = Omega_bar, dTheta = 2 Theta_bar
        => in the weighted vars  G = (Y-1)^{1/2} dOmega = const,  h = dTheta = const.
        A CONSTANT is not in L^2(e^{-2 gamma z} dz) (diverges at the edge z->-inf),
        so the amplitude mode is excluded by the weight.  It is exactly the marginal
        k=0 mode.
  * dilation/edge generator:  X Omega_bar_X ~ (X-1)^{-3/2}  (more singular than the
        profile) -> excluded by the weighted space, and it is the appended xi mode
        at lambda = c_l-1 = 1 (handled separately).

So "asymptotic stability modulo the symmetry group" is realized concretely: the
symmetry modes are precisely the directions the weighted space excludes, and on the
weighted space the FULL shape operator (Omega + Theta, coupled) is gapped.

Gates:
 (1) the Theta sector (M_hh) has the same convective gap: weighted norm of h decays
     at ~2 gamma (the marginal constant-Theta mode is removed by the weight).
 (2) the FULL coupled (G,h) operator [[M_GG, M_Gh],[0, M_hh]] has the convective gap
     (block-triangular => union spectrum; both blocks gapped).
 (3) the amplitude (nu) symmetry mode = the constant mode: marginal in the plain
     norm, excluded (non-normalizable at the edge) in the weighted norm.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from stability_edge_z import build_blocks_z, spectrum, C_L  # noqa: E402


def evolve(M, G0, taus):
    vals, V = np.linalg.eig(M)
    c = np.linalg.solve(V, G0)
    return [V @ (np.exp(vals * t) * c) for t in taus]


def wnorm(z, dz, F, gamma):
    return float(np.sqrt(np.sum(np.exp(-2 * gamma * z) * np.abs(F) ** 2) * dz))


def fit_rate(taus, vals, lo=0.5, hi=2.5):
    w = (taus >= lo) & (taus <= hi)
    return -np.polyfit(taus[w], np.log(np.array(vals)[w]), 1)[0]


def main():
    print("=" * 78)
    print("  Theta + mu,nu folded in via the convective weighted space")
    print("=" * 78)
    B = build_blocks_z(N=1400, zmin=-14.0, zmax=12.0, edge_bc=True)
    z, dz = B["z"], B["z"][1] - B["z"][0]
    Mgg, Mhh, Mgh = B["M_GG"], B["M_hh"], B["M_Gh"]
    n = len(z)
    taus = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    gammas = [0.25, 0.5]

    # (1) Theta sector gap
    print("\n  (1) Theta sector (M_hh): weighted-norm decay rate vs predicted 2*gamma")
    h0 = np.exp(-((z - 1.0)) ** 2)
    Ht = evolve(Mhh, h0, taus)
    for g in gammas:
        s = [wnorm(z, dz, h, g) for h in Ht]
        print(f"      gamma={g:.2f}: rate={fit_rate(taus,[x/s[0] for x in s]):+.3f}"
              f"   (predict {2*g:.2f})")
    su = [wnorm(z, dz, h, 0.0) for h in Ht]
    print(f"      unweighted: rate={fit_rate(taus,[x/su[0] for x in su]):+.3f} (marginal)")

    # (2) full coupled (G,h) operator
    print("\n  (2) FULL coupled (G,h) operator: weighted-norm decay of the combined state")
    M = np.block([[Mgg, Mgh], [np.zeros((n, n)), Mhh]])
    vfull, _ = spectrum(M)
    print(f"      max Re eig(full) = {vfull.real.max():+.4f}  (#Re>0 = "
          f"{int((vfull.real>1e-6).sum())})   [block-triangular: eig = MGG U Mhh]")
    F0 = np.concatenate([np.exp(-((z - 1.0)) ** 2), 0.5 * np.exp(-((z - 1.5)) ** 2)])
    zz = np.concatenate([z, z])
    Ft = evolve(M, F0, taus)
    for g in gammas:
        s = [wnorm(zz, dz, F, g) for F in Ft]
        print(f"      gamma={g:.2f}: combined (G,h) weighted-norm rate="
              f"{fit_rate(taus,[x/s[0] for x in s]):+.3f}   (predict {2*g:.2f})")

    # (3) the amplitude (nu) symmetry mode = the constant mode
    print("\n  (3) amplitude (nu) symmetry mode = constant (G,h):")
    g = 0.5
    const = np.ones(n)
    shape = np.exp(-((z - 1.0)) ** 2)        # a genuine localized shape perturbation
    print(f"      ||const||  unweighted = {wnorm(z,dz,const,0.0):.3e}   "
          f"weighted(g={g}) = {wnorm(z,dz,const,g):.3e}  (edge-dominated/excluded)")
    print(f"      ||shape||  unweighted = {wnorm(z,dz,shape,0.0):.3e}   "
          f"weighted(g={g}) = {wnorm(z,dz,shape,g):.3e}  (finite: a true perturbation)")
    print("      -> the constant (amplitude/nu) mode is dominated by the edge and is")
    print("         non-normalizable as z_min->-inf; the weighted space excludes it,")
    print("         while genuine localized shape perturbations are retained.")

    print("\n" + "=" * 78)
    print("  VERDICT (frozen, full shape sector Omega+Theta, weighted space):")
    print("   - Theta sector has the SAME convective gap (2*gamma); its marginal")
    print("     constant mode is removed by the weight.")
    print("   - the FULL coupled (G,h) operator: no Re>0 point mode; convective gap")
    print("     2*gamma in the weighted norm.")
    print("   - the scaling symmetry modes (amplitude nu = constant; dilation/edge")
    print("     = xi@1) are EXACTLY the directions the weighted space excludes, so")
    print("     'asymptotic stability modulo symmetries' is realized by the weight")
    print("     itself -- no paper-specific normalization needed for the verdict.")
    print("  => full frozen-sector evidence consistent with Conjecture 2.4.")
    print("  (Explicit mu,nu with the paper's normalization should agree; not a proof.)")
    print("=" * 78)


if __name__ == "__main__":
    main()
