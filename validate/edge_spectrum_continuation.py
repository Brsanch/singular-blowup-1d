#!/usr/bin/env python3
"""Spectral-side continuation for the frozen Omega-shape generator M_GG (Conjecture 2.4):
what is the rightmost spectrum of M_GG, and why a direct eigenvalue computation cannot
certify the asymptotic-stability GAP. (docs/desk/edge_layer_estimate.md, ledger OBS-04.)

ROBUST findings (reproduced below; cross-checked over z_min, N, far-field weight beta, and
an edge g-mass damping kappa):

  (A) max Re eig(M_GG) ~ -1e-3 everywhere, and -> 0 as N -> inf.  NO eigenvalue with Re > 0
      appears at any z_min in [-18,-12], N in [1400,2800].  But the maximum is MARGINAL: it
      is the transport ESSENTIAL spectrum sitting on Re=0, not a strict spectral gap.

  (B) That near-0 cluster is the ESSENTIAL (transport) spectrum, NOT the inadmissible
      (Y-1)^{-1/2} edge branch: an edge g-mass damping -kappa*diag(e^{-z}) -- which crushes
      the +40 field-of-values abscissa to negative at kappa=1e-3 (OBS-03) -- does NOT move
      the near-0 eigenvalues even at kappa=1e-2.  They are delocalized and come in conjugate
      pairs.  So the field-of-values +40 (kappa-crushable non-normality on the inadmissible
      branch) and the essential spectrum on Re=0 (kappa-stable, genuine) are DISTINCT objects.

  (C) Consequence: the asymptotic-stability gap (Re spec <= -delta) is irreducibly a
      WEIGHTED-space property.  The transport essential spectrum is on Re=0 unweighted; only a
      weight shifts it left, and that weight is exactly what diverges at the edge (cap_scaffold
      + docs/desk/edge_layer_estimate.md).  A direct eigenvalue computation of the unweighted
      M_GG therefore CANNOT exhibit the gap -- it shows the marginal essential spectrum at Re=0.

NOT CLAIMED: a clean value for "the rightmost genuine point eigenvalue".  A naive smoothness
(hi-k) filter does NOT robustly separate the discrete point spectrum from the edge/essential
modes on this strongly non-normal operator (the selection flips with weighted-vs-unweighted
scoring).  And the discrete Hilbert mis-resolves the profile at the edge (OBS-03), so the
edge-decay / c1=H(g)(1) readouts are unreliable there.  This script characterizes the spectrum;
it does NOT numerically prove stability or B=0 -- that is the analytic leap.
Pure numpy eig -> no apfsd/panic risk.
"""
import sys
sys.path.insert(0, "/Volumes/4TB SD/ClaudeCode/singular-blowup-1d/src")
import numpy as np
from stability_edge_z import build_blocks_z, spectrum


print("=" * 90)
print("(A) max Re eig(M_GG) vs (z_min, N)  --  Re>0 count, and the value -> 0 as N->inf")
print("    zmin    N=1400      N=2000      N=2800     (#Re>0 in parens)")
for zmin in (-12., -14., -16., -18.):
    row = []
    for N in (1400, 2000, 2800):
        v, _ = spectrum(build_blocks_z(N=N, zmin=zmin, zmax=12., beta=0.0, edge_bc=True)["M_GG"])
        npos = int((v.real > 1e-6).sum())
        row.append(f"{v.real.max():+.5f}({npos})")
    print(f"    {zmin:+.0f}    " + "   ".join(row))

print("\n(B) near-0 cluster under edge g-mass damping -kappa*diag(e^{-z})")
print("    (kappa-stable => delocalized essential spectrum, NOT the inadmissible edge branch)")
for (zmin, N) in [(-12., 2000), (-16., 2000)]:
    B = build_blocks_z(N=N, zmin=zmin, zmax=12., beta=0.0, edge_bc=True)
    pen = np.diag(np.exp(-B["z"]))
    print(f"    zmin={zmin:+.0f} N={N}: rightmost 4 Re eig(M_GG):")
    for kappa in (0.0, 1e-4, 1e-2):
        v, _ = spectrum(B["M_GG"] - kappa * pen)
        print(f"      kappa={kappa:.0e}:  " + "  ".join(f"{x.real:+.4f}" for x in v[:4]))

print("\n" + "=" * 90)
print("READ: no Re>0 mode anywhere, but max Re is the MARGINAL essential spectrum on Re=0")
print("(kappa-stable, -> 0 with N) -- the repo's 'benign -1e-3'.  The stability GAP needs the")
print("weight; the weight breaks at the edge => the open analytic leap (not a numerical question).")
print("=" * 90)
