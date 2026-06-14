#!/usr/bin/env python3
"""Independent reproduction of the edge-layer workflow's decisive claims:
(1) the full energy abscissa A_full collapses as A_full * exp(zmin/2) = const
    (continuum profile-amplitude (Y-1)^{-1/2} rate, NOT a grid rate);
(2) dropping a fixed-OFFSET edge layer leaves a NEGATIVE bulk abscissa, depth-independent;
(3) an infinitesimal g-mass penalty (kappa) collapses the +40 (it lives on the
    g-unbounded branch).
Uses the repo operator directly. Pure numpy eigh -> no apfsd/panic risk.
"""
import sys, os
sys.path.insert(0, "/Volumes/4TB SD/ClaudeCode/singular-blowup-1d/src")
import numpy as np
from stability_edge_z import build_blocks_z

def phi_prime(z, gamma, z0, s):
    return 2.0 * gamma * 0.5 * (1.0 + np.tanh((z - z0) / s))

def energy_op(B, gamma=0.4, z0=-6.0, s=2.0):
    z, dz, Qop = B["z"], B["dz"], B["Q_op"]
    pp = phi_prime(z, gamma, z0, s)
    phi = np.cumsum(pp) * dz
    sw = np.exp(-0.5 * phi)
    Ft = (sw[:, None] * (0.5 * Qop)) / sw[None, :]
    S = -np.diag(pp) + 0.5 * (Ft + Ft.T)
    return S, z, sw, pp

print("=" * 78)
print("(1) SELF-SIMILAR COLLAPSE  A_full * exp(zmin/2) = const ?  (N=1400)")
print("    zmin    A_full      A_full*exp(zmin/2)    bulk(z>-6)   bulk(z>-4)")
for zmin in (-10., -12., -14., -16.):
    B = build_blocks_z(N=1400, zmin=zmin, zmax=12.0, beta=0.0, edge_bc=True)
    S, z, sw, pp = energy_op(B)
    A = float(np.linalg.eigvalsh(S)[-1])
    coll = A * np.exp(zmin / 2)
    b6 = float(np.linalg.eigvalsh(S[np.ix_(z > -6, z > -6)])[-1])
    b4 = float(np.linalg.eigvalsh(S[np.ix_(z > -4, z > -4)])[-1])
    print(f"    {zmin:+.0f}   {A:+9.3f}      {coll:+.4f}            {b6:+.4f}      {b4:+.4f}")

print("\n(2) N-STABILITY of the collapse constant at fixed zmin=-12")
print("    N       A_full      A_full*exp(-6)")
for N in (900, 1400, 2000, 2800):
    B = build_blocks_z(N=N, zmin=-12.0, zmax=12.0, beta=0.0, edge_bc=True)
    S, z, sw, pp = energy_op(B)
    A = float(np.linalg.eigvalsh(S)[-1])
    print(f"    {N:<6d}  {A:+9.3f}      {A*np.exp(-6.):+.4f}")

print("\n(3) g-MASS PENALTY: add kappa * diag(e^{-z}) (penalize g=e^{-z/2}G unboundedness)")
print("    the +40 lives on the g-unbounded branch -> small kappa should crush it")
B = build_blocks_z(N=1400, zmin=-12.0, zmax=12.0, beta=0.0, edge_bc=True)
S, z, sw, pp = energy_op(B)
A0 = float(np.linalg.eigvalsh(S)[-1])
print(f"    kappa=0:      A_full={A0:+.3f}")
for kappa in (1e-4, 1e-3, 1e-2):
    Sk = S - kappa * np.diag(np.exp(-z))   # subtract: penalize edge (e^{-z} huge at edge)
    Ak = float(np.linalg.eigvalsh(Sk)[-1])
    print(f"    kappa={kappa:.0e}:  A_full={Ak:+.3f}   ({100*(A0-Ak)/A0:.0f}% reduction)")

print("\n(4) DISCRETE HILBERT mis-resolves the profile at the deepest node?")
print("    H((Y-1)^{-1/2}) should be 0 in continuum (G-R 3.222.2); check deepest node")
from stability_edge_z import hilbert_support_z
for zmin in (-10., -12., -14., -16.):
    Nn = 1400
    z = np.linspace(zmin, 12.0, Nn); Y = 1.0 + np.exp(z)
    Ymax = Y[-1] + 0.5*(Y[-1]-Y[-2])
    Hz = hilbert_support_z(Y, Ymax)
    g = (Y - 1.0) ** -0.5
    Hg = Hz @ g
    print(f"    zmin={zmin:+.0f}: H(profile)[node0]={Hg[0]:+.1f}   exp(-zmin/2)={np.exp(-zmin/2):.1f}   ratio={Hg[0]/np.exp(-zmin/2):+.3f}")
print("=" * 78)
