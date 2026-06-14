#!/usr/bin/env python3
"""Linchpin numerics for the edge Mellin/Wiener-Hopf analysis (docs/desk/edge_mellin_symbol.md).

In z=log(Y-1) the half-line (Y>1) Hilbert transform is a CONVOLUTION with pv-kernel
1/(e^{z-z'}-1), hence a Fourier multiplier / Mellin symbol acting on powers (Y-1)^p:

  (1)  H[(Y-1)^p] = cot(pi p) * (Y-1)^p      (Mellin symbol H(p)=cot(pi p), zero at p=-1/2
                                              reproduces the profile identity H[(Y-1)^{-1/2}]=0)
  (2)  (1/2)Q_op[(Y-1)^zeta] = -[tan(pi zeta)/(2 zeta+1)] (Y-1)^{zeta-1/2}
       (the -1/2 exponent shift = the profile amplitude; coefficient = 0 at integer zeta
        => regular edge modes lambda=-2 zeta<=0; pole at zeta=1/2 => forced-constant +40 resonance)

Valid in the convergence strip (-1/2 < zeta < 1/2 for the forcing; -1<p<0 for H); outside it the
single-power readout exits the strip (far-field domination) and the symbol must be analytically
continued -- which is the Wiener-Hopf content, NOT a discretization error. Pure numpy -> no panic.
"""
import sys
sys.path.insert(0, "/Volumes/4TB SD/ClaudeCode/singular-blowup-1d/src")
import numpy as np
from stability_edge_z import hilbert_support_z, build_blocks_z

print("=" * 84)
print("(1) Hilbert Mellin symbol:  H[(Y-1)^p]/(Y-1)^p  =?  cot(pi p)   (strip -1<p<0)")
z = np.linspace(-16.0, 16.0, 6000); Y = 1.0 + np.exp(z)
Ymax = Y[-1] + 0.5 * (Y[-1] - Y[-2]); H = hilbert_support_z(Y, Ymax); ym1 = Y - 1.0
win = (z > -6) & (z < 4)
print("   p       cot(pi p)     measured      match")
for p in (-0.8, -0.7, -0.5, -0.3, -0.2):
    r = float(np.median((H @ ym1 ** p) / ym1 ** p)[win] if False else np.median(((H @ ym1 ** p) / ym1 ** p)[win]))
    t = 1.0 / np.tan(np.pi * p)
    print(f"   {p:+.2f}    {t:+.5f}      {r:+.5f}     {'YES' if abs(r-t)<0.05*(abs(t)+0.2) else 'edge-of-strip'}")

print("\n(2) Forcing symbol:  (1/2)Q_op[(Y-1)^zeta] =?  -tan(pi z)/(2z+1) * (Y-1)^{z-1/2}")
B = build_blocks_z(N=6000, zmin=-16.0, zmax=16.0, beta=0.0, edge_bc=False)
z2, Qop, ym = B["z"], B["Q_op"], B["Y"] - 1.0; win2 = (z2 > -6) & (z2 < 2)
print("   zeta   out-power   z-0.5    coef        pred           note")
for zeta in (-0.3, 0.0, 0.3, 0.5, 0.7):
    FG = 0.5 * (Qop @ ym ** zeta)
    powfit = float(np.polyfit(z2[win2], np.log(np.abs(FG[win2]) + 1e-300), 1)[0])
    coef = float(np.median(FG[win2] / (ym[win2] ** (zeta - 0.5))))
    d = 2 * zeta + 1
    pred = -np.tan(np.pi * zeta) / d if abs(d) > 1e-9 else np.inf
    note = ("in-strip" if -0.5 < zeta < 0.5 else ("RESONANCE" if abs(zeta-0.5) < 1e-6 else "out-of-strip (continue)"))
    ps = f"{pred:+.4f}" if np.isfinite(pred) and abs(pred) < 1e6 else "  pole/inf"
    print(f"   {zeta:+.2f}   {powfit:+.3f}     {zeta-0.5:+.2f}    {coef:+.4f}     {ps}      {note}")
print("=" * 84)
