#!/usr/bin/env python3
"""
GATE: convective-weight structure pack for M_GG = -2 d_z + (1/2)Q_op  (LEM-06 attack,
docs/desk/convective_weight_spectrum.md, 2026-07-15).

Machine-verifies, with NO discrete Hilbert transform anywhere (continuum quadratures and
Gamma/Mellin identities only):

  S1  z-form functional equation:  (lam+2s) X(s) = c(s+1/2) X(s+1/2),
      X(s) = Gamma(2s+lam)Gamma(2s+2)Gamma(s)Gamma(1-s),  c(s) = -tan(pi s)/(2s+1).
  S2  ladder symbol by continuum quadrature:  K e^{sz} = c(s) e^{(s-1/2)z}.
  S3  defect sign:  K(bump) -> -Ghat(1/2)/(2 pi) at the edge  (domain constraint of L2_sigma).
  S4  residue:  c(s+1/2) ~ 1/(2 pi s)  at s = 0.
  S5  closed form of the ladder coefficients:  Pi_n(s) = X(s)/[X(s+n/2)(lam+2s+n)].
  S6  rank-one defect identity:  M G_p = lam G_p - Gamma(1+lam)*1   (G_p = X-inversion), plus
      far-field constant Gamma(lam) and e^z edge decay.
  S7  ladder resolvent solves (lam - M) G = F  (defect-closed; F = e^{-cosh(z-z0)}), with
      Ghat(1/2) = 0 and vanishing far-field constant (residue<->value identity).
  S8  energy no-go mechanism: exact diagonalization w = e^{-z/4}u gives multiplier
      Re c(sigma+1/4+it); Re c(sigma+1/4) > 0 for sigma in (0,1/2) (positive band), while
      Gaussian bumps have Re-form < 0 with exact e^{-a/2} depth scaling.
  S9  vertical-line Barnes rotation obstruction: integrand log-slope -> pi/4 (the reason the
      Section-8 resummation contour must be bent).

Exit 0 iff all PASS.
"""
import sys

import mpmath as mp
import numpy as np
from scipy.special import gamma as G_

FAILURES = []


def report(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}  {detail}")
    if not ok:
        FAILURES.append(name)


SIGMA = 0.4

# ---------------------------------------------------------------- mp helpers
mp.mp.dps = 30


def c_mp(s):
    return -mp.tan(mp.pi * s) / (2 * s + 1)


def X_mp(s, lam):
    return mp.gamma(2 * s + lam) * mp.gamma(2 * s + 2) * mp.gamma(s) * mp.gamma(1 - s)


# ================================================================ S1
print("S1: functional equation of the 4-Gamma solution")
worst = mp.mpf(0)
for lam in [mp.mpf('0.3'), mp.mpc('0.5', '1.5'), mp.mpf('-0.37'), mp.mpc('2.3', '-0.8')]:
    for s in [mp.mpf('0.2'), mp.mpc('0.35', '2.0'), mp.mpc('0.1', '-5.0')]:
        lhs = (lam + 2 * s) * X_mp(s, lam)
        rhs = c_mp(s + mp.mpf('0.5')) * X_mp(s + mp.mpf('0.5'), lam)
        worst = max(worst, abs(lhs - rhs) / (abs(lhs) + abs(rhs)))
report("S1 functional equation", worst < mp.mpf('1e-24'), f"worst rel {mp.nstr(worst, 3)}")

# ================================================================ S2
print("S2: ladder symbol by continuum quadrature (pv Hilbert + causal)")


def H_power(p, Y):
    a = Y - 1
    eps = mp.mpf('1e-12')
    f = lambda u: u ** p / (a - u)
    i1 = mp.quad(f, [0, a * (1 - eps)])
    i2 = mp.quad(f, [a * (1 + eps), 2 * a])
    U = a * mp.mpf('1e8')
    i3 = mp.quad(f, [2 * a, 256 * a, mp.mpf('1e6') * a, U])
    i3 += U ** p / p + a * U ** (p - 1) / (p - 1)
    return (i1 + i2 + i3) / mp.pi


worst = mp.mpf(0)
for s in [mp.mpf('0.30'), mp.mpf('0.45')]:
    p = s - mp.mpf('0.5')
    z0 = mp.mpf('0.7')
    val = mp.quad(lambda zp: mp.e ** zp * H_power(p, 1 + mp.e ** zp), [mp.mpf('-60'), z0])
    KG = mp.mpf('0.5') * mp.e ** (-z0) * val
    pred = c_mp(s) * mp.e ** ((s - mp.mpf('0.5')) * z0)
    worst = max(worst, abs(KG - pred) / abs(pred))
report("S2 ladder symbol", worst < mp.mpf('1e-8'), f"worst rel {mp.nstr(worst, 3)}")

# ================================================================ S3
print("S3: defect sign on a bump (edge constant of K G)")
z_c = mp.mpf('1.0')
Gb = lambda z: mp.e ** (-(z - z_c) ** 2)
Ghat_half = mp.quad(lambda z: Gb(z) * mp.e ** (-z / 2), [-mp.inf, mp.inf])


def H_g_bump(zp):
    ezp = mp.e ** zp
    f = lambda zpp: mp.e ** (-zpp / 2) * Gb(zpp) * mp.e ** zpp / (ezp - mp.e ** zpp)
    lo, hi = mp.mpf('-9'), mp.mpf('11')
    eps = mp.mpf('1e-9')
    if zp < lo or zp > hi:
        return mp.quad(f, [lo, hi]) / mp.pi
    return (mp.quad(f, [lo, zp - eps]) + mp.quad(f, [zp + eps, hi])) / mp.pi


kg_edge = mp.mpf('0.5') * mp.e ** 9 * mp.quad(lambda zp: mp.e ** zp * H_g_bump(zp),
                                              [mp.mpf('-45'), mp.mpf('-9')])
pred = -Ghat_half / (2 * mp.pi)
report("S3 defect sign", abs(kg_edge - pred) / abs(pred) < mp.mpf('1e-3'),
       f"K(bump)(-9) = {mp.nstr(kg_edge, 6)} vs -Ghat(1/2)/2pi = {mp.nstr(pred, 6)}")

# ================================================================ S4
print("S4: residue of c(s+1/2) at s=0")
eps = mp.mpf('1e-8')
res = eps * c_mp(eps + mp.mpf('0.5'))
report("S4 residue 1/(2pi)", abs(res - 1 / (2 * mp.pi)) < mp.mpf('1e-6'),
       f"{mp.nstr(res, 8)} vs {mp.nstr(1/(2*mp.pi), 8)}")

# ================================================================ S5
print("S5: closed form Pi_n(s) = X(s)/[X(s+n/2)(lam+2s+n)]")
lam = mp.mpf('0.3')
s = mp.mpf(str(SIGMA)) + 1j * mp.mpf('0.42')
Pi = 1 / (lam + 2 * s)
worst = mp.mpf(0)
for n in range(10):
    closed = X_mp(s, lam) / (X_mp(s + mp.mpf(n) / 2, lam) * (lam + 2 * s + n))
    worst = max(worst, abs(Pi - closed) / abs(closed))
    Pi = Pi * c_mp(s + mp.mpf(n + 1) / 2) / (lam + 2 * s + (n + 1))
report("S5 closed form", worst < mp.mpf('1e-22'), f"worst rel {mp.nstr(worst, 3)}")

# ================================================================ S6/S7 (numpy hybrid)
print("S6/S7: defect identity + ladder resolvent (vectorized continuum operator)")
TAU_MAX, DTAU = 14.0, 0.02
ZLO, ZHI, DZ = -45.0, 45.0, 0.006
taus = np.arange(-TAU_MAX, TAU_MAX + DTAU / 2, DTAU)
S_line = SIGMA + 1j * taus
zg = np.arange(ZLO, ZHI + DZ / 2, DZ)
ezg = np.exp(zg)
trapz = getattr(np, "trapezoid", None) or np.trapz


def c_np(s):
    return -np.tan(np.pi * s) / (2 * s + 1)


def X_np(s, lam):
    return G_(2 * s + lam) * G_(2 * s + 2) * G_(s) * G_(1 - s)


def invert(samples, z):
    out = np.empty(len(z), dtype=complex)
    B = 4000
    for i0 in range(0, len(z), B):
        zc = z[i0:i0 + B]
        ph = np.exp((SIGMA + 1j * taus)[None, :] * zc[:, None])
        out[i0:i0 + B] = trapz(samples[None, :] * ph, taus, axis=1) / (2 * np.pi)
    return out


def K_apply(Gz):
    gz = np.exp(-zg / 2) * Gz
    num = gz * ezg
    with np.errstate(divide='ignore', invalid='ignore'):
        pv_hi = ZHI - np.log(np.abs(ezg - np.exp(ZHI)))
        pv_lo = ZLO - np.log(np.abs(ezg - np.exp(ZLO)))
    H = np.empty_like(Gz, dtype=complex)
    for i in range(len(zg)):
        a = ezg[i]
        d = a - ezg
        d[i] = 1.0
        sub = (num - num[i]) / d
        sub[i] = 0.0
        I_sub = DZ * (np.sum(sub) - 0.5 * (sub[0] + sub[-1]))
        H[i] = (I_sub + num[i] * (pv_hi[i] - pv_lo[i]) / a) / np.pi
    H[0] = H[1]
    H[-1] = H[-2]
    integ = ezg * H
    cum = np.concatenate([[0.0], np.cumsum((integ[1:] + integ[:-1]) / 2 * DZ)])
    return 0.5 * np.exp(-zg) * cum


lamf = 0.35
Xs = X_np(S_line, lamf)
Gp = invert(Xs, zg)
MGp = -2 * invert(S_line * Xs, zg) + K_apply(Gp)
pred = lamf * Gp - G_(1 + lamf)
mask = (zg > -35) & (zg < 20)
err = np.max(np.abs(MGp[mask] - pred[mask])) / (np.max(np.abs(pred[mask])) + abs(G_(1 + lamf)))
report("S6 defect identity", err < 3e-3, f"rel resid {err:.2e}")
i_far = np.argmin(abs(zg - 40))
i30 = np.argmin(abs(zg + 30))
i28 = np.argmin(abs(zg + 28))
ok_far = abs(Gp[i_far] - G_(lamf)) < 0.02 * abs(G_(lamf))
ok_edge = abs(abs(Gp[i30] / Gp[i28]) - np.exp(-2)) < 1e-3
report("S6 far-const/edge-decay", ok_far and ok_edge,
       f"G_p(40)={Gp[i_far]:.4f} vs Gamma={G_(lamf):.4f}; edge ratio {abs(Gp[i30]/Gp[i28]):.5f}")

lamf, z0 = 0.3, -6.0
vg = np.arange(-16, 16, 0.002)
ECOSH = np.exp(-np.cosh(vg))


def Fhat_batch(svals):
    out = np.empty(len(svals), dtype=complex)
    B = 200
    for i0 in range(0, len(svals), B):
        sc = svals[i0:i0 + B]
        mat = ECOSH[None, :] * np.exp(-sc[:, None] * vg[None, :])
        out[i0:i0 + B] = trapz(mat, vg, axis=1)
    return out * np.exp(-svals * z0)


def G0hat_line(svals, nmax=80):
    total = np.zeros(len(svals), dtype=complex)
    Pi = 1.0 / (lamf + 2 * svals)
    for n in range(nmax):
        term = Pi * Fhat_batch(svals + n / 2.0)
        total += term
        if np.max(np.abs(term)) < 1e-17 * np.max(np.abs(total) + 1e-300):
            break
        Pi = Pi * c_np(svals + (n + 1) / 2.0) / (lamf + 2 * svals + (n + 1))
    return total


G0s = G0hat_line(S_line)
# G0hat(1/2) via the exact closed form Pi_n(1/2) = X(1/2)/[X(1/2+n/2)(lam+1+n)]
# (odd n vanish; even n finite). A finite +-delta limit is analytically wrong at
# depth (dps-stable garbage — see validate/edge_convective_resolvent_scan.py).
G0half = 0.0 + 0.0j
Xhalf = X_np(np.array([0.5 + 0j]), lamf)[0]
for n in range(0, 60, 2):
    sn = np.array([0.5 + n / 2.0 + 0j])
    term = Xhalf / (X_np(sn, lamf)[0] * (lamf + 1 + n)) * complex(Fhat_batch(sn)[0])
    G0half += term
    if n > 6 and abs(term) < 1e-17 * (abs(G0half) + abs(term)):
        break
fix = G0half / (2 * np.pi * G_(1 + lamf))
Gs = G0s - fix * X_np(S_line, lamf)

s_chk = np.array([SIGMA + 0.7345j])
g_s = complex((G0hat_line(s_chk) - fix * X_np(s_chk, lamf))[0])
g_sh = complex((G0hat_line(s_chk + 0.5) - fix * X_np(s_chk + 0.5, lamf))[0])
fh = complex(Fhat_batch(s_chk)[0])
line_resid = abs((lamf + 2 * s_chk[0]) * g_s - c_np(s_chk[0] + 0.5) * g_sh - fh) / abs(fh)
report("S7 line equation", line_resid < 1e-12, f"rel resid {line_resid:.2e}")

G = invert(Gs, zg)
MG = -2 * invert(S_line * Gs, zg) + K_apply(G)
F = np.exp(-np.cosh(np.clip(zg - z0, -700, 700)))
pde = np.max(np.abs((lamf * G - MG - F)[mask])) / np.max(F)
report("S7 resolvent PDE", pde < 5e-3, f"rel resid {pde:.2e}")
report("S7 no far-field constant", abs(G[i_far]) < 1e-3,
       f"G(40) = {abs(G[i_far]):.2e} (residue<->value identity)")

# ================================================================ S8
print("S8: energy no-go — positive band of Re c(sigma+1/4+it), dissipative bumps")
band_ok = all(complex(c_mp(mp.mpf(x))).real > 0 for x in ['0.52', '0.65', '0.74'])
report("S8 positive band Re c(sigma+1/4) > 0", band_ok,
       f"Re c(0.52,0.65,0.74) = {[float(mp.re(c_mp(mp.mpf(x)))) for x in ['0.52','0.65','0.74']]}")

mp.mp.dps = 20
phi = lambda v: mp.e ** (-v ** 2)


def K_form(a):
    f = lambda t: c_mp(mp.mpf(str(SIGMA)) + 1j * t) * (mp.sqrt(mp.pi) * mp.e ** (-t ** 2 / 4)) \
        * mp.e ** (-(mp.mpf(str(SIGMA)) + 1j * t) * a + mp.mpf(str(SIGMA)) * a)
    def Ku(z):
        return mp.quad(lambda t: f(t) * mp.e ** ((mp.mpf(str(SIGMA)) + 1j * t - mp.mpf('0.5')) * z),
                       [-30, -6, 0, 6, 30]) / (2 * mp.pi)
    num = mp.quad(lambda z: Ku(z) * (mp.e ** (mp.mpf(str(SIGMA)) * z) * phi(z - a)) * mp.e ** (-2 * mp.mpf(str(SIGMA)) * z),
                  [a - 8, a, a + 8])
    return num


f3, f5 = K_form(mp.mpf('-3')), K_form(mp.mpf('-5'))
scale_ok = abs(f3 * mp.e ** mp.mpf('-1.5') - f5 * mp.e ** mp.mpf('-2.5')) / abs(f3 * mp.e ** mp.mpf('-1.5')) < mp.mpf('1e-5')
neg_ok = mp.re(f3) < 0
report("S8 bump form e^{-a/2}-scaling + dissipative sign", bool(scale_ok and neg_ok),
       f"form*e^{{a/2}}: {mp.nstr(f3*mp.e**mp.mpf('-1.5'), 7)} vs {mp.nstr(f5*mp.e**mp.mpf('-2.5'), 7)}")

# ================================================================ S9
print("S9: vertical-line Barnes obstruction (log-slope -> pi/4)")
mp.mp.dps = 30
lam = mp.mpf('0.3')
s = mp.mpf(str(SIGMA)) + 1j * mp.mpf('0.42')
theta = mp.mpf('0.3')
SH = {}


def shape(sv):
    key = mp.nstr(sv, 20)
    if key not in SH:
        SH[key] = mp.quad(lambda v: mp.e ** (-mp.cosh(v) - sv * v), [-16, -8, -4, 0, 4, 8, 16])
    return SH[key]


def integrand_mag(t):
    nu = -theta + 1j * t
    return abs((mp.pi / mp.sin(mp.pi * nu)) * shape(s + nu / 2) * mp.e ** (-(s + nu / 2) * mp.mpf('-6'))
               / (X_mp(s + nu / 2, lam) * (lam + 2 * s + nu)))


slope = mp.log(integrand_mag(mp.mpf('34')) / integrand_mag(mp.mpf('30'))) / 4
report("S9 growth slope in (0.6, pi/4+0.06)", bool(mp.mpf('0.6') < slope < mp.pi / 4 + mp.mpf('0.06')),
       f"slope {mp.nstr(slope, 4)} (pi/4 = {mp.nstr(mp.pi/4, 4)})")

# ================================================================
print()
if FAILURES:
    print(f"GATE FAIL: {FAILURES}")
    sys.exit(1)
print("GATE PASS: all convective-weight structure checks verified")
