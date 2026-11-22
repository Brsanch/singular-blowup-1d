"""Task 5/7: a physical-time check of (1.1) that genuinely depends on V_G.
V_G(tau) is computed from the *physical* profile by the formula in (1.1); Q by the real-space kernel on the
physical profile (no dilation-covariance shortcut); G_tau by 4th-order differences along the s-grid.
The horizon ODE dr/ds = -1 + mut(V)/2 is integrated with the true V and with mutated V."""
import math, sys
import numpy as np
from lib import *
from solver import *
from scipy.interpolate import CubicSpline
A = float(sys.argv[1]) if len(sys.argv) > 1 else 0.05
sc = G0_maxnorm_scale()
S, nst = 1.2, 2400
s = Solver(-64.0, 64.0, 4096, sponge=10.0)
z = s.g.z; dz = s.g.dz
s.setup(S / nst)
h = A * h0_family(s) / sc
Gs = [s.N(h)]
for k in range(nst):
    h = s.step(h); Gs.append(s.N(h))
Gs = np.array(Gs)
Aq = Q_real_matrix(z)

def VG_s(G, r):  # s-clock form used for the ODE
    arg = r * r * np.exp(-z)
    ln_term = np.where(arg > 1e6, 2.0 * math.log(r) - z + np.log1p(1.0 / np.maximum(arg, 1e-300)), np.log1p(arg))
    return float(-(1.0 / (math.pi * r)) * np.sum(np.exp(0.5 * z) * ln_term * G) * dz)

def V_phys(Gp):  # formula of (1.1) on the physical profile, y = e^zeta
    return float(-(1.0 / math.pi) * np.sum(np.exp(0.5 * z) * np.log1p(np.exp(-z)) * Gp) * dz)

print(f"A={A}:  residual of (1.1) with V from the physical profile; columns: s, r, tau, V_phys, V_ode, max|res| over v in [e^-3,e^3], max|V v G_v|, max|NL terms|")
for lab, mut in (("true V_G", lambda v: v), ("ODE uses 3 V_G", lambda v: 3 * v), ("ODE uses -V_G", lambda v: -v), ("ODE uses 0.05 sin(7r)", None)):
    r = np.zeros(nst + 1); r[0] = 1.0; Vo = np.zeros(nst + 1)
    f = (lambda G, rr: 0.05 * math.sin(7 * rr)) if mut is None else (lambda G, rr: mut(VG_s(G, rr)))
    kend = nst
    for k in range(nst):
        Vo[k] = f(Gs[k], r[k]); k1 = -1 + 0.5 * Vo[k]
        rp = r[k] + s.ds * k1
        if rp <= 1e-3: kend = k; break
        k2 = -1 + 0.5 * f(Gs[k + 1], rp)
        r[k + 1] = r[k] + 0.5 * s.ds * (k1 + k2)
        if r[k + 1] <= 1e-3: kend = k; break
    tau = np.concatenate([[0], np.cumsum(0.5 * s.ds * (1 / r[:kend] + 1 / r[1:kend + 1]))])
    for s0 in (0.3, 0.6, 0.9):
        k0 = int(round(s0 / s.ds))
        if k0 + 2 > kend: continue
        zc = z[(z > -3) & (z < 3)][::8]
        def Gphys(k, zeta):
            return CubicSpline(z, Gs[k])(zeta + 2 * math.log(r[k]))
        vals = [Gphys(k0 + j, zc) for j in (-2, -1, 1, 2)]
        dGds = (vals[0] - 8 * vals[1] + 8 * vals[2] - vals[3]) / (12 * s.ds)
        G_tau = r[k0] * dGds
        Gp_full = CubicSpline(z, Gs[k0])(np.clip(z + 2 * math.log(r[k0]), z[0], z[-1]))
        Gp_full[z + 2 * math.log(r[k0]) < z[0] + 1] = 0.0
        Qp = Q_real(Gp_full, z, Aq)
        DGp = s.D(Gp_full)
        sel = np.isin(z, zc)
        Vp = V_phys(Gp_full)
        res = G_tau + (2 - Vp + Qp[sel]) * DGp[sel] - 0.5 * Qp[sel] * (1 + Gp_full[sel])
        print(f"  {lab:22s} s={s0:.1f} r={r[k0]:.4f} tau={tau[k0]:.3f}  V_phys={Vp:+.5e} V_ode={Vo[k0]:+.5e}  max|res|={np.max(np.abs(res)):.2e}"
              f"   |V v G_v|={abs(Vp)*np.max(np.abs(DGp[sel])):.2e}  |Q v G_v|,|Q G/2|={np.max(np.abs(Qp[sel]*DGp[sel])):.2e},{np.max(np.abs(0.5*Qp[sel]*Gp_full[sel])):.2e}")
