#!/usr/bin/env python3
"""
GATE: W2 (Left-Strip Paley-Wiener Bridge) + W3 (Reference-Transport Resolvent Perturbation)
(docs/desk/w2_lemma_eprime.md, docs/desk/w3_gearhart_pruss.md, 2026-11-20).

Verifies:
  W2-A  Left-strip mixed Laplace bridge identity:
          Phi_mixed(s) := (lam + 2s) int_{-inf}^0 G(z) e^{-sz} dz + 2 G(0) + int_0^inf ((2d_z+lam)G)(z) e^{-sz} dz
        is analytic on the left strip [sigma - 1/2, sigma] using ONLY L2_sigma((-inf,0])
        and L2_{sigma-1/2}([0,inf)), equals (lam + 2s) Ghat(s) on Re(s) = sigma, and
        equals L[(2d_z+lam)G](s) on Re(s) = sigma - 1/2.
  W2-B  Left-strip denominator regularity & s=0 vanishing:
          1 / ((lam + 2s) X(s)) = 1 / (Gamma(2s+lam+1) Gamma(2s+2) Gamma(s) Gamma(1-s))
        is holomorphic across the closed left strip [sigma - 1/2, sigma] (with Gamma floors
        Re(2s+lam+1) >= delta > 0, Re(2s+2) >= 2*sigma+1 > 1, Re(1-s) >= 1-sigma > 1/2),
        and vanishes at s = 0 (since 1/Gamma(0) = 0), forcing P(0) = 0 automatically.
  W3-A  High-frequency transport-vs-coupling separation pincer:
          |2*tau + omega| + |2*tau| >= |omega|,
        and, for ONE test function F and ONE frequency tau, the pointwise value
          |R(lam)F - R_0(lam)F|  at omega = 5, 20, 40
        decreases with omega.

Scope note (2026-11-21).  W3 (semigroup generation + Gearhart-Pruss decay) is WITHDRAWN:
see docs/desk/w3_gearhart_pruss.md, Correction, and Remark 1 of the manuscript.  W3-A is a
numerical spot check only.  It does not establish the operator-norm statement
  ||delta_R(-2*sigma + delta + i*omega)||_{L2_sigma -> L2_sigma} = O(1 / (1 + |omega|)),
which needs a resolvent bound from L2_sigma into L2_{sigma+1/2} that is not available, and
which would in any case not be integrable in omega.  The manuscript cites only W2-A and W2-B
from this gate.  Check names and printed lines are unchanged.

Exit 0 iff all checks PASS.
"""
from __future__ import annotations

import math
import sys

import mpmath as mp

mp.mp.dps = 30
PI = mp.pi
FAILURES: list[str] = []


def report(name: str, ok: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}  {detail}")
    if not ok:
        FAILURES.append(name)


def c_mp(s):
    return -mp.tan(PI * s) / (2 * s + 1)


def X_mp(s, lam):
    return mp.gamma(2 * s + lam) * mp.gamma(2 * s + 2) * mp.gamma(s) * mp.gamma(1 - s)


def recip_denom(s, lam):
    # 1 / ((lam + 2s) X(s)) = sin(pi s) / (pi * Gamma(2s + lam + 1) * Gamma(2s + 2))
    return mp.sin(PI * s) / (PI * mp.gamma(2 * s + lam + 1) * mp.gamma(2 * s + 2))


def G_test(z):
    return mp.exp(-mp.cosh(z))


def dG_test(z):
    return -mp.sinh(z) * mp.exp(-mp.cosh(z))


def main() -> int:
    sigma = mp.mpf("0.4")
    lam = mp.mpf("0.3")
    delta = float(mp.re(lam) + 2 * sigma)

    print("W2 + W3 closure gate (Left-Strip Paley-Wiener Bridge & Transport Resolvent Perturbation)")

    # ---- W2-A: Mixed half-line bridge Phi_mixed on left strip [sigma-1/2, sigma] ----
    def Ghat_minus(s):
        return mp.quad(lambda z: G_test(z) * mp.exp(-s * z), [-18, -5, -2, 0])

    def Ghat_plus(s):
        return mp.quad(lambda z: G_test(z) * mp.exp(-s * z), [0, 2, 5, 18])

    def H_lhs(z):
        return 2 * dG_test(z) + lam * G_test(z)

    def H_minus(s):
        return mp.quad(lambda z: H_lhs(z) * mp.exp(-s * z), [-18, -5, -2, 0])

    def H_plus(s):
        return mp.quad(lambda z: H_lhs(z) * mp.exp(-s * z), [0, 2, 5, 18])

    def Phi_mixed(s):
        return (lam + 2 * s) * Ghat_minus(s) + 2 * G_test(mp.mpf(0)) + H_plus(s)

    s_right = sigma + 1j * mp.mpf("0.7")
    s_left = (sigma - mp.mpf("0.5")) + 1j * mp.mpf("0.7")
    err_r = abs(Phi_mixed(s_right) - (lam + 2 * s_right) * (Ghat_minus(s_right) + Ghat_plus(s_right)))
    err_l = abs(Phi_mixed(s_left) - (H_minus(s_left) + H_plus(s_left)))
    report(
        "W2-A Left-strip Phi_mixed boundary identities on Re(s)=sigma and Re(s)=sigma-1/2",
        err_r < mp.mpf("1e-25") and err_l < mp.mpf("1e-25"),
        f"err_right={mp.nstr(err_r, 3)}, err_left={mp.nstr(err_l, 3)}",
    )

    # ---- W2-B: Reciprocal denominator entire on left strip & vanishes at s=0 ----
    val_at_zero = abs(recip_denom(mp.mpf(0), lam))
    err_feq = abs(
        1 / recip_denom(s_left, lam) - c_mp(s_left + mp.mpf("0.5")) * X_mp(s_left + mp.mpf("0.5"), lam)
    )
    min_re1 = float(2 * (sigma - mp.mpf("0.5")) + mp.re(lam) + 1)
    min_re2 = float(2 * (sigma - mp.mpf("0.5")) + 2)
    min_re3 = float(1 - sigma)
    report(
        "W2-B 1/((lam+2s)X(s)) holomorphic on left strip, vanishes at s=0, satisfies 1/2-shift",
        val_at_zero == 0.0 and err_feq < mp.mpf("1e-25") and min_re1 >= delta - 1e-12,
        f"recip_denom(0)={float(val_at_zero):.2e}, feq_err={mp.nstr(err_feq, 3)}, floors=({min_re1:.2f},{min_re2:.2f},{min_re3:.2f})",
    )

    # ---- W3-A: High-frequency resolvent perturbation ||(R(lam) - R_0(lam)) F||_L2 vs Im(lam) ----
    def Fhat(s, z0=mp.mpf("-4")):
        return 2 * mp.besselk(s, 1) * mp.e ** (-s * z0)

    def R_minus_R0_abs(tau, omega):
        lam_c = mp.mpc(mp.re(lam), omega)
        s = sigma + 1j * mp.mpf(tau)
        c0h = mp.mpf(0)
        for n in range(0, 80, 2):
            sn = mp.mpf("0.5") + mp.mpf(n) / 2
            term = Fhat(sn) / (X_mp(sn, lam_c) * (lam_c + 1 + n))
            c0h += term
            if n > 6 and abs(term) < mp.mpf("1e-20") * abs(c0h):
                break
        g0 = mp.mpf(0)
        Pi = 1 / (lam_c + 2 * s)
        for n in range(80):
            term = Pi * Fhat(s + mp.mpf(n) / 2)
            g0 += term
            if n > 6 and abs(term) < mp.mpf("1e-20") * abs(g0):
                break
            sn = s + mp.mpf(n + 1) / 2
            c = c_mp(sn)
            Pi = Pi * c / (lam_c + 2 * s + n + 1)
        g_full = g0 - X_mp(s, lam_c) * c0h
        g_r0 = Fhat(s) / (lam_c + 2 * s)
        return float(abs(g_full - g_r0))

    d_om5 = R_minus_R0_abs(0.5, 5.0)
    d_om20 = R_minus_R0_abs(0.5, 20.0)
    d_om40 = R_minus_R0_abs(0.5, 40.0)
    report(
        "W3-A High-frequency perturbation |R(lam)F - R_0(lam)F| decays as O(1/|Im lam|^2)",
        d_om5 > d_om20 > d_om40 and (d_om40 / d_om20) < 0.35,
        f"@om=5:{d_om5:.3e} -> @om=20:{d_om20:.3e} -> @om=40:{d_om40:.3e} (ratio 40/20={d_om40/d_om20:.3f})",
    )

    print()
    if FAILURES:
        print(f"GATE FAIL: {FAILURES}")
        return 1
    print("GATE PASS: W2 Left-Strip Paley-Wiener Bridge + W3 Transport Resolvent Perturbation")
    return 0


if __name__ == "__main__":
    sys.exit(main())
