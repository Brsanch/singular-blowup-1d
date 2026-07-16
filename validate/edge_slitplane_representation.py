#!/usr/bin/env python3
"""
GATE: slit-plane Plemelj representation of the ladder resolvent + the dressed-Schur
closure pieces (docs/desk/convective_weight_spectrum.md Section 10, 2026-07-16).

P1  Representation: in w = e^{4 pi i s} the coboundary C(s) - C(s+1/2) = Phi(s) is an
    additive jump problem across the ray arg w = 4 pi sigma; the SPLIT-regularized
    Cauchy transform (plain kernel on the w'->0 wing, one subtraction on the infinity
    wing — Muskhelishvili class-h bookkeeping forced by |Phi| ~ e^{5 pi |tau'|/2}) gives

        Ghat(s) = X(s) [ D(e^{4 pi i s}) - D(1) ]      (constant-freedom-free)

    == the down-ladder defect-closed solution, at strip-interior points, to ~1e-20.
    Jump orientation: D_- - D_+ = Phi in the tau'-parametrization (sign -1 below).
    This formula consumes ONLY line values of Fhat => it defines the solution for
    arbitrary L2_sigma data (the minus-Hardy wall has an explicit classical formula).

P2  Four-quadrant Schur envelope: the dressed split kernel
        |X(tau)|/|X(tau')| * ker(tau,tau') / |lam + 2 s'|
    (ker = plain 4pi|w'|/|w'-w| on tau'>0, subtracted 4pi|w|/|w'-w| on tau'<0)
    satisfies  dressed * e^{+pi|tau-tau'|}  bounded-and-decaying in EVERY quadrant
    (the 4pi - 3pi = pi margin; the former "corner failure" was the unsubtracted
    kernel used on the wrong wing).

P3  High-frequency corner benign (OBS-08 pin): for F = e^{-cosh(z-z0)} e^{i Omega z}
    the dressing gain/transport -> 1 as |Omega| grows (coupling c(sigma+i Omega) ~
    1/(2|Omega|)); pins dressing(z0=-6, Om=-15) ~ 0.9968.

Exit 0 iff all PASS.
"""
import sys

import mpmath as mp

mp.mp.dps = 40
SIGMA = mp.mpf('0.4')
LAM = mp.mpf('0.3')
Z0 = mp.mpf('-6')
FAILURES = []


def report(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}  {detail}")
    if not ok:
        FAILURES.append(name)


def c_sym(s):
    return -mp.tan(mp.pi * s) / (2 * s + 1)


def X(s, lam=LAM):
    return mp.gamma(2*s + lam) * mp.gamma(2*s + 2) * mp.gamma(s) * mp.gamma(1 - s)


SHAPE = {}


def K1(nu):
    key = mp.nstr(nu, 30)
    if key not in SHAPE:
        SHAPE[key] = 2 * mp.besselk(nu, 1)
    return SHAPE[key]


def Fhat(s, Om=mp.mpf(0)):
    return K1(s - 1j*Om) * mp.e**(-(s - 1j*Om) * Z0)


# ---------------------------------------------------------------- P1
def Phi(taup):
    s = SIGMA + 1j * taup
    return Fhat(s) / (X(s) * (LAM + 2*s))


def wray(taup):
    return mp.e**(4j * mp.pi * SIGMA) * mp.e**(-4 * mp.pi * taup)


def D(w, sign=-1):
    T = mp.mpf(12)
    f_plus = lambda tp: Phi(tp) * (-4 * mp.pi * wray(tp)) / (wray(tp) - w)
    f_minus = lambda tp: Phi(tp) * (-4 * mp.pi * wray(tp)) * (1/(wray(tp) - w) - 1/wray(tp))
    return sign * (mp.quad(f_plus, [0, 1, 3, 6, T]) + mp.quad(f_minus, [-T, -6, -3, -1, 0])) / (2j * mp.pi)


def G0hat(s, Om=mp.mpf(0), nmax=300):
    tol = mp.mpf(10) ** (-30)
    total = mp.mpf(0)
    Pi = 1 / (LAM + 2*s)
    for n in range(nmax):
        term = Pi * Fhat(s + mp.mpf(n)/2, Om)
        total += term
        if n > 6 and abs(term) < tol * (abs(total) + abs(term)):
            break
        Pi = Pi * c_sym(s + mp.mpf(n+1)/2) / (LAM + 2*s + (n+1))
    return total


def C0half(Om=mp.mpf(0), nmax=400):
    half = mp.mpf('0.5')
    total = mp.mpf(0)
    for n in range(0, nmax, 2):
        sn = half + mp.mpf(n)/2
        term = Fhat(sn, Om) / (X(sn) * (LAM + 1 + n))
        total += term
        if n > 6 and abs(term) < mp.mpf(10)**(-30) * (abs(total) + abs(term)):
            break
    return total


print("P1: slit-plane Cauchy representation == ladder solution (strip interior)")
Dhalf = D(mp.mpf(1))
c0h = C0half()
worst = mp.mpf(0)
for s_test in [SIGMA + mp.mpf('0.06'), SIGMA + mp.mpf('0.10') + 1j*mp.mpf('0.9'),
               SIGMA + mp.mpf('0.45')]:
    lad = X(s_test) * (G0hat(s_test)/X(s_test) - c0h)
    cau = X(s_test) * (D(mp.e**(4j * mp.pi * s_test)) - Dhalf)
    worst = max(worst, abs(lad - cau) / (abs(lad) + abs(cau)))
report("P1 representation", worst < mp.mpf('1e-15'), f"worst rel {mp.nstr(worst, 3)}")

# ---------------------------------------------------------------- P2
print("P2: four-quadrant dressed-Schur envelope (x e^{pi|tau-tau'|} bounded, decaying)")


def dressed(t, tp):
    w, wp = wray(t), wray(tp)
    ker = 4 * mp.pi * abs(wp) / abs(wp - w)
    if tp < 0:
        ker *= abs(w) / abs(wp)
    return abs(X(SIGMA + 1j*t)) / abs(X(SIGMA + 1j*tp)) * ker / abs(LAM + 2*(SIGMA + 1j*tp))


ok2 = True
for name, pts in [("corner tp<t<0", [(-1, -3), (-4, -12)]),
                  ("tp>t>0", [(1, 3), (4, 12)]),
                  ("0<tp<t", [(3, 1), (12, 4)]),
                  ("t<tp<0", [(-3, -1), (-12, -4)]),
                  ("opposite", [(2, -2), (-4, 4)])]:
    vals = [dressed(mp.mpf(t), mp.mpf(tp)) * mp.e**(mp.pi * abs(mp.mpf(t) - mp.mpf(tp)))
            for t, tp in pts]
    ok_q = all(v < 1 for v in vals) and vals[-1] < vals[0] + mp.mpf('1e-10')
    ok2 &= bool(ok_q)
    print(f"    {name:>14}: {'  '.join(mp.nstr(v, 3) for v in vals)}  {'ok' if ok_q else 'FAIL'}")
report("P2 Schur envelope", ok2)

# ---------------------------------------------------------------- P3
print("P3: high-frequency corner benign (dressing -> 1)")
mp.mp.dps = 60          # cancellation here is only ~1e7; clear the dps-40 cache
SHAPE.clear()
Om = mp.mpf('-15')


def gain_ratio():
    # uniform Simpson tau in [Om-12, 12]
    lo, hi = float(Om) - 12, 12.0
    m = int((hi - lo) / 0.125)
    if m % 2 == 1:
        m += 1
    h = (hi - lo) / m
    taus = [mp.mpf(lo) + j*h for j in range(m+1)]
    ws = [h/3 * (1 if j in (0, m) else (4 if j % 2 == 1 else 2)) for j in range(m+1)]
    c0 = C0half(Om, nmax=800)
    nG2 = mp.mpf(0)
    nA2 = mp.mpf(0)
    for t, w in zip(taus, ws):
        s = SIGMA + 1j*t
        g = G0hat(s, Om, nmax=500) - c0 * X(s)
        fh = Fhat(s, Om)
        nG2 += w * abs(g)**2
        nA2 += w * abs(fh / (LAM + 2*s))**2
    return mp.sqrt(nG2 / nA2)


dr = gain_ratio()
report("P3 dressing(z0=-6, Om=-15) in (0.98, 1.005)", bool(mp.mpf('0.98') < dr < mp.mpf('1.005')),
       f"dressing = {mp.nstr(dr, 5)}")

print()
if FAILURES:
    print(f"GATE FAIL: {FAILURES}")
    sys.exit(1)
print("GATE PASS: slit-plane representation + Schur-closure pieces verified")
