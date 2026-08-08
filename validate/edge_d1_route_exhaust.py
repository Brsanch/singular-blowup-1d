#!/usr/bin/env python3
"""
W1-D1 route exhaust probe (2026-08-07 goal).

Does NOT claim D1 proved, LEM-06 closed, or Conjecture 2.4.
Does NOT use Schur S_off or OBS-07 gains as a D1 discharge.

What this measures / certifies numerically:

  S1  Difference identity (load-bearing for G):
        C0(s) - C0(1/2)  ==  D(w(s)) - D(1)
      on strip-interior test points (rel err << 1e-12).
      Absolute equality D(1) == C0(1/2) is FALSE (rel O(0.1)).

  S2  Resolvent G matches from either difference formula at a generic s
      (not the removable s=1/2 point).

  S3  Slit kernel for the scalar D(1) as a functional of F on the line:
        K(tau) = (-1/(2 j pi)) * kappa_split(tau, w=1)
                 / (X(sigma+i tau) * (lam + 2(sigma+i tau)))
      |K| e^{+pi |tau|} is poly-bounded (decaying in practice);
      truncated ||K||_L2 saturates (K in L2 as NUMERIC-FORECAST).

  S4  |D(1)| <= ||K||_2 ||F||_2 on the bilateral bump family
      (CS sanity; not an operator-norm proof).

Exit 0 iff S1–S3 structural checks pass. S4 is reported, not gated.
False-friend ban: no Schur/OBS-07 discharge lines.
"""
from __future__ import annotations

import sys

import mpmath as mp

mp.mp.dps = 35
SIGMA = mp.mpf("0.4")
LAM = mp.mpf("0.3")
Z0 = mp.mpf("-6")
PI = mp.pi
FAILURES: list[str] = []


def report(name: str, ok: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}  {detail}")
    if not ok:
        FAILURES.append(name)


def X(s, lam=LAM):
    return (
        mp.gamma(2 * s + lam)
        * mp.gamma(2 * s + 2)
        * mp.gamma(s)
        * mp.gamma(1 - s)
    )


_SHAPE: dict = {}


def K1(nu):
    key = mp.nstr(nu, 28)
    if key not in _SHAPE:
        _SHAPE[key] = 2 * mp.besselk(nu, 1)
    return _SHAPE[key]


def Fhat(s, z0=Z0):
    return K1(s) * mp.e ** (-s * z0)


def wray(tp):
    return mp.e ** (4j * PI * SIGMA) * mp.e ** (-4 * PI * tp)


def Phi(tp, z0=Z0):
    s = SIGMA + 1j * tp
    return Fhat(s, z0) / (X(s) * (LAM + 2 * s))


def D_of(w, z0=Z0, T=mp.mpf(12)):
    """Split-regularized Cauchy transform at w (gate formula, sign=-1)."""

    def f_plus(tp):
        return Phi(tp, z0) * (-4 * PI * wray(tp)) / (wray(tp) - w)

    def f_minus(tp):
        return Phi(tp, z0) * (-4 * PI * wray(tp)) * (
            1 / (wray(tp) - w) - 1 / wray(tp)
        )

    return (
        (-1)
        * (
            mp.quad(f_plus, [0, 1, 3, 6, T])
            + mp.quad(f_minus, [-T, -6, -3, -1, 0])
        )
        / (2j * PI)
    )


def G0(s, z0=Z0, nmax=250):
    total = mp.mpf(0)
    Pi = 1 / (LAM + 2 * s)
    for n in range(nmax):
        term = Pi * Fhat(s + mp.mpf(n) / 2, z0)
        total += term
        if n > 6 and abs(term) < mp.mpf("1e-22") * (abs(total) + abs(term)):
            break
        sn = s + mp.mpf(n + 1) / 2
        c = -mp.tan(mp.pi * sn) / (2 * sn + 1)
        Pi = Pi * c / (LAM + 2 * s + (n + 1))
    return total


def C0_half(z0=Z0, nmax=300):
    """Even-n closed form at s=1/2 (safe; never +-delta)."""
    total = mp.mpf(0)
    for n in range(0, nmax, 2):
        sn = mp.mpf("0.5") + mp.mpf(n) / 2
        term = Fhat(sn, z0) / (X(sn) * (LAM + 1 + n))
        total += term
        if n > 6 and abs(term) < mp.mpf("1e-22") * (abs(total) + abs(term)):
            break
    return total


def kappa_split(tp, w=mp.mpf(1)):
    wp = wray(tp)
    if tp >= 0:
        return (-4 * PI * wp) / (wp - w)
    return (-4 * PI * wp) * (1 / (wp - w) - 1 / wp)


def K_slit(tp, lam=LAM):
    """Kernel: D(1) = int K(tp) F(sigma+i tp) d tp  (for the split definition)."""
    s = SIGMA + 1j * tp
    return (
        mp.mpf(-1)
        / (2j * PI)
        * kappa_split(tp)
        / (X(s, lam) * (lam + 2 * s))
    )


def K_L2(T=mp.mpf(15), n=1200, lam=LAM):
    h = mp.mpf(2 * T) / n
    acc = mp.mpf(0)
    for j in range(n + 1):
        tp = -T + j * h
        w = mp.mpf("0.5") if j in (0, n) else mp.mpf(1)
        acc += w * abs(K_slit(tp, lam)) ** 2
    return mp.sqrt(acc * h)


def F_L2(z0, T=mp.mpf(15), n=800):
    h = mp.mpf(2 * T) / n
    acc = mp.mpf(0)
    for j in range(n + 1):
        tp = -T + j * h
        w = mp.mpf("0.5") if j in (0, n) else mp.mpf(1)
        s = SIGMA + 1j * tp
        acc += w * abs(Fhat(s, z0)) ** 2
    return mp.sqrt(acc * h)


# ---------------------------------------------------------------------------
print("W1-D1 route exhaust probe (no D1/LEM-06 proof claim)")
print(f"  sigma={SIGMA}  lam={LAM}  z0={Z0}  dps={mp.mp.dps}")
print()

# S1 difference identity
print("S1: difference identity C0(s)-C0(1/2) == D(w)-D(1); absolute D-id false")
D1 = D_of(mp.mpf(1))
c0h = C0_half()
abs_rel = abs(D1 - c0h) / (abs(D1) + abs(c0h))
report(
    "S1a absolute D(1)!=C0(1/2) (rel>1e-2)",
    abs_rel > mp.mpf("1e-2"),
    f"rel={mp.nstr(abs_rel, 4)}  D1={mp.nstr(D1, 6)}  C0={mp.nstr(c0h, 6)}",
)

worst_diff = mp.mpf(0)
const_spread = mp.mpf(0)
base_const = c0h - D1
for s in [
    SIGMA + mp.mpf("0.06"),
    SIGMA + mp.mpf("0.12") + 1j * mp.mpf("0.7"),
    SIGMA + mp.mpf("0.35"),
]:
    w = mp.e ** (4j * PI * s)
    C0s = G0(s) / X(s)
    dC = C0s - c0h
    dD = D_of(w) - D1
    rel = abs(dC - dD) / (abs(dC) + abs(dD) + mp.mpf("1e-30"))
    worst_diff = max(worst_diff, rel)
    const_spread = max(const_spread, abs((C0s - D_of(w)) - base_const))
report(
    "S1b difference identity",
    worst_diff < mp.mpf("1e-12"),
    f"worst rel={mp.nstr(worst_diff, 3)}  const_spread={mp.nstr(const_spread, 3)}",
)

# S2 G match at generic s (avoid s with w=1)
print("S2: G from ladder difference == G from slit difference")
s_g = SIGMA + mp.mpf("0.18") + 1j * mp.mpf("0.4")
G_lad = X(s_g) * (G0(s_g) / X(s_g) - c0h)
G_cau = X(s_g) * (D_of(mp.e ** (4j * PI * s_g)) - D1)
grel = abs(G_lad - G_cau) / (abs(G_lad) + abs(G_cau) + mp.mpf("1e-30"))
report("S2 G match", grel < mp.mpf("1e-12"), f"rel={mp.nstr(grel, 3)}")

# S3 kernel decay / L2
print("S3: slit kernel K for D(1) — decay margin and ||K||_2 saturation")
scaled_ok = True
for t in [2, 4, 6, 8, 10]:
    tp = mp.mpf(t)
    scaled = abs(K_slit(tp)) * mp.e ** (PI * abs(tp))
    # must not explode; allow poly growth but we see decay
    if scaled > mp.mpf(10):
        scaled_ok = False
    print(f"    tau={t:+d}  |K|e^{{pi|t|}}={mp.nstr(scaled, 5)}")
    tp = -tp
    scaled = abs(K_slit(tp)) * mp.e ** (PI * abs(tp))
    if scaled > mp.mpf(10):
        scaled_ok = False
    print(f"    tau={float(tp):+g}  |K|e^{{pi|t|}}={mp.nstr(scaled, 5)}")

# Mass of |K|^2 is near tau=0; tails e^{-2 pi |tau|} — compare mid windows.
# Very large |tau| underflows |X|^{-1} at modest dps and can bias T-huge norms.
k5 = K_L2(T=mp.mpf(5), n=500)
k10 = K_L2(T=mp.mpf(10), n=900)
k14 = K_L2(T=mp.mpf(14), n=1100)
sat = abs(k14 - k10) / (k14 + mp.mpf("1e-30"))
report(
    "S3a |K|e^{pi|t|} bounded (forecast margin)",
    scaled_ok,
    f"||K||_5={mp.nstr(k5, 6)}  ||K||_10={mp.nstr(k10, 6)}  ||K||_14={mp.nstr(k14, 6)}",
)
report(
    "S3b ||K||_2 saturates (T=10 vs 14)",
    sat < mp.mpf("5e-3"),
    f"rel change={mp.nstr(sat, 4)}",
)

# S4 CS sanity (report only)
print("S4: CS sanity |D(1)| <= ||K||_2 ||F||_2 on bumps (report only)")
Kn = k14
print(f"    ||K||_2 ≈ {mp.nstr(Kn, 7)}")
for z0 in [mp.mpf(0), mp.mpf(-4), mp.mpf(2)]:
    d1 = abs(D_of(mp.mpf(1), z0=z0))
    fn = F_L2(z0)
    ratio = d1 / (Kn * fn + mp.mpf("1e-30"))
    print(
        f"    z0={mp.nstr(z0, 3)}  |D1|={mp.nstr(d1, 5)}  "
        f"||F||={mp.nstr(fn, 5)}  |D1|/(||K||||F||)={mp.nstr(ratio, 4)}"
    )

print()
print("DISCIPLINE: Schur S_off / OBS-07 gains are NOT used as D1 discharge.")
print("CLAIM BAN: this probe does not prove D1, LEM-06, or PDE blowup.")
print()
if FAILURES:
    print(f"PROBE FAIL: {FAILURES}")
    sys.exit(1)
print("PROBE PASS: structural S1–S3 (difference identity + K forecast in L2)")
sys.exit(0)
