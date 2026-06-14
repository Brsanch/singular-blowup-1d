#!/usr/bin/env python3
"""Argument-principle eigenvalue COUNT of M_GG in the right half-plane (Conjecture 2.4).

Two robust counts, cross-checked, of #{point eigenvalues of M_GG with Re lambda >= eps}:
  (A) direct: eig(M_GG, edge_bc=True), count Re >= eps  (nested eps; vs N, z_min; + g-penalty).
  (B) contour winding: N_count = (1/2pi i) oint tr((M - lambda I)^{-1}) dlambda  around a box
      {eps <= Re lambda <= Lambda, |Im lambda| <= T}.  Counts genuine eigenvalues inside the box,
      robust to clustering; the essential spectrum sits on Re~0 (LEFT of the box) so is excluded.
Expected: 0 (no unstable Omega-shape point mode), robust => numerical evidence (with the
discrete-Hilbert edge-bias caveat for large Re lambda) for the Wiener-Hopf index = 0.
Pure numpy/scipy -> no apfsd/panic risk.
"""
import sys
sys.path.insert(0, "/Volumes/4TB SD/ClaudeCode/singular-blowup-1d/src")
import numpy as np
from stability_edge_z import build_blocks_z

# ---------------- (A) direct eigenvalue count ----------------
print("=" * 88)
print("(A) #eig(M_GG, edge_bc=True) with Re lambda >= eps   (nested eps; vs N, z_min; + g-penalty)")
print("    zmin  N     #Re>=.05  #Re>=.25  #Re>=.5  #Re>=1   maxRe     | +g-penalty(1e-2) #Re>=.05")
for zmin in (-12., -14., -16.):
    for N in (1400, 2000):
        B = build_blocks_z(N=N, zmin=zmin, zmax=12., beta=0.0, edge_bc=True)
        v = np.linalg.eigvals(B["M_GG"]).real
        c = [int((v >= e).sum()) for e in (0.05, 0.25, 0.5, 1.0)]
        # g-penalty kills the inadmissible (Y-1)^{-1/2} edge branch
        vg = np.linalg.eigvals(B["M_GG"] - 1e-2 * np.diag(np.exp(-B["z"]))).real
        cg = int((vg >= 0.05).sum())
        print(f"   {zmin:+.0f}  {N:<5d} {c[0]:>6d}   {c[1]:>7d}  {c[2]:>7d}  {c[3]:>6d}   {v.max():+.4f}   |   {cg:>6d}")

# ---------------- (B) contour-integral winding count ----------------
print("\n(B) contour count  N = (1/2pi i) oint tr((M-lambda)^{-1}) dlambda  around the box")
print("    {eps<=Re<=Lambda, |Im|<=T}  (genuine eigenvalues enclosed; essential spectrum on Re~0 excluded)")
def contour_count(M, eps=0.1, Lam=6.0, T=8.0, npts=400):
    n = M.shape[0]; I = np.eye(n)
    # rectangle corners (ccw): (eps,-T)->(Lam,-T)->(Lam,T)->(eps,T)->(eps,-T)
    corners = [eps - 1j*T, Lam - 1j*T, Lam + 1j*T, eps + 1j*T, eps - 1j*T]
    total = 0.0 + 0j
    per = max(40, npts // 4)
    for a, b in zip(corners[:-1], corners[1:]):
        ts = (np.arange(per) + 0.5) / per
        lam = a + (b - a) * ts
        dl = (b - a) / per
        for L in lam:
            # tr((M-L I)^{-1}) via solve
            X = np.linalg.solve(M - L * I, I)
            total += np.trace(X) * dl
    return total / (2j * np.pi)

for (zmin, N) in [(-12., 1400), (-14., 1400)]:
    B = build_blocks_z(N=N, zmin=zmin, zmax=12., beta=0.0, edge_bc=True)
    val = contour_count(B["M_GG"], eps=0.1, Lam=6.0, T=8.0, npts=480)
    print(f"   zmin={zmin:+.0f} N={N}: contour count = {val.real:+.3f} (Im {val.imag:+.1e})  -> round = {int(round(val.real))}")

print("=" * 88)
print("READ: all #Re>=eps == 0 and contour count ~ 0 => NO unstable Omega-shape point mode,")
print("robust in eps, N, z_min, and under the inadmissible-branch g-penalty. (Continuum out-of-strip")
print("Re lambda>=1 still needs the analytic WH index; the discrete Hilbert biases the deep edge.)")
