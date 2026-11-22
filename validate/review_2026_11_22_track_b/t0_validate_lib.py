"""Validate reviewer's toolkit before using it."""
import math, sys
import numpy as np
import mpmath as mp
from lib import *

# (1) N(m) vs mpmath
mp.mp.dps = 30
for m in (0.25 + 3j, 0.625 - 7.5j, -0.25 + 0.3j, 0.625 + 40j, 0.25 - 0.0j):
    ref = mp.cos(mp.pi * m) * mp.gamma(2 * m + 2) / mp.pi
    got = Nsym(np.array([m]))[0]
    print("N", m, abs(got - complex(ref)) / abs(complex(ref)))
print("N'(1/2) =", (Nsym(np.array([0.5 + 1e-6]))[0] - Nsym(np.array([0.5 - 1e-6]))[0]).real / 2e-6)

# (2) real-space Q vs closed form for G = w/(b+w)^2  (a != 0)
b = 2.0
for (zmin, zmax, n) in ((-40, 60, 2000), (-40, 60, 4000), (-60, 80, 4000)):
    g = Grid(zmin, zmax, n)
    z, w = g.z, g.w
    G0 = w / (b + w) ** 2
    qex = np.where(w < 1e-8, -1 / (2 * b**1.5) + 3 * w / (4 * b**2.5),
                   (1 / (2 * math.sqrt(b) * w)) * np.log1p(w / b) - 1 / (math.sqrt(b) * (b + w)))
    q = Q_real(G0, z)
    for lo, hi in ((-30, -15), (-15, 0), (0, 15), (15, 30)):
        mk = (z > lo) & (z < hi)
        print(f"Qreal grid[{zmin},{zmax}] n={n} region({lo},{hi}): abs err {np.max(np.abs(q-qex)[mk]):.2e}  rel-to-max|q| {np.max(np.abs(q-qex)[mk])/np.max(np.abs(qex[mk])):.2e}")

# (2b) q0 exact itself vs mpmath PV definition at a few w (independent of both codes)
def q_def(wv):
    f = lambda y: mp.sqrt(y) / (b + y) ** 2 * mp.log(abs(wv / y - 1))
    return mp.quad(f, [0, wv / 2, wv, 2 * wv, 50 * wv + 50, mp.inf]) / (mp.pi * wv)
for wv in (0.01, 1.0, 30.0):
    ex = (1 / (2 * math.sqrt(b) * wv)) * math.log1p(wv / b) - 1 / (math.sqrt(b) * (b + wv))
    print("q def vs closed form", wv, float(q_def(mp.mpf(wv))), ex)

# (3) conjugation (1/2) Q (N h) = -N_r (u h) with REAL-space Q, h = u^3 e^{-u}  (a = 0)
g = Grid(-48, 48, 4096)
ops = conj_ops(g)
h0 = g.u**3 * np.exp(-g.u)
G = ops["N"](h0)
qc = ops["Q"](G)
qr = Q_real(G, g.z)
for lo, hi in ((-30, -15), (-15, 0), (0, 15), (15, 30)):
    mk = (g.z > lo) & (g.z < hi)
    print(f"conj vs real Q region({lo},{hi}): {np.max(np.abs(qc-qr)[mk]):.2e}  (max|Q| {np.max(np.abs(qr[mk])):.2e})")
print("a(G) quad:", a_quad(G, g.z), " max|G|", np.max(np.abs(G)))
# (4) closed form check: N phi_* = (1/2) u^2 e^{-u}  (derived by reviewer: Ghat = Gamma(2m+2))
g = Grid(-200, 120, 2**15)
ops = conj_ops(g, sig=0.0)
phi = g.u / (g.u**2 + 1)
F = ops["N"](phi)
mk = (g.z > -20) & (g.z < 20)
print("N phi_* vs u^2 e^{-u}/2 (line 0):", np.max(np.abs(F - 0.5 * g.u**2 * np.exp(-g.u))[mk]))
