#!/usr/bin/env python3
"""
W1-D1 dual probe — empirical defect functional + formal Hardy dual kernel.

Does NOT claim D1 proved. Does NOT re-run Schur/OBS-07 as a discharge.

What this measures (NUMERIC-FORECAST only):

  (P1) On the dense class F_z0(s) = 2 K_s(1) e^{-s z0}
       (Laplace of e^{-cosh(z-z0)}; entire in s), the even-ladder
           D1(F) = sum_{k>=0} a_k F(1/2+k)
           a_k = 1 / [X(1/2+k) (lam + 1 + 2k)]
       vs ||F||_L2(line Re s=sigma). Ratio R(z0,lam) = |D1|/||F|| is a
       lower bound on the op-norm of any continuous extension.

  (P2) Formal H^2-Cauchy dual on the line:
           psi(tau) = (1/2pi) sum_k a_k / (sigma + i tau - s_k)
       If D1 extended continuously through H^2 of Re s > sigma, this would
       be the Riesz representer on that subspace. Diagnostics:
         - leading large-|tau| coeff C_inf := sum_k a_k  (psi ~ C_inf/(2pi i tau))
         - truncated ||psi||_L2 on [-T,T] growth in T  (1/tau => log; 1/tau^2 => bound)

  (P3) delta-scaling: R at bulk z0 vs delta = Re(lam)+2*sigma.

Exit 0 always after printing a structured report (this is a probe, not a gate).
Pass/fail lines are diagnostic only.

Ref: docs/desk/w1_plemelj_stirling_writeup.md §3.5 (W1-D1),
     docs/desk/convective_weight_spectrum.md Prop F / §8.
"""
from __future__ import annotations

import math
import sys

import mpmath as mp

DPS = 40
mp.mp.dps = DPS
SIGMA = mp.mpf("0.4")
PI = mp.pi


def X(s, lam):
    return (
        mp.gamma(2 * s + lam)
        * mp.gamma(2 * s + 2)
        * mp.gamma(s)
        * mp.gamma(1 - s)
    )


def a_k(k: int, lam):
    """Even-ladder coefficient at s = 1/2 + k (n = 2k)."""
    s = mp.mpf("0.5") + k
    n = 2 * k
    return 1 / (X(s, lam) * (lam + 1 + n))


def F_bump(s, z0):
    """Entire transform of e^{-cosh(z-z0)}."""
    return 2 * mp.besselk(s, 1) * mp.e ** (-s * z0)


def D1_even(lam, z0, kmax: int = 40):
    total = mp.mpf(0)
    for k in range(kmax):
        s = mp.mpf("0.5") + k
        term = a_k(k, lam) * F_bump(s, z0)
        total += term
        if k > 4 and abs(term) < mp.mpf("1e-30") * (abs(total) + 1):
            break
    return total


def F_L2_line(lam_unused, z0, sigma=SIGMA, tau_max=20, n_tau=400):
    """||F||_L2(d tau / 2pi)^{1/2} on Re s = sigma via trapezoid."""
    # Plancherel convention matching desk: measure d tau / 2pi
    h = mp.mpf(2 * tau_max) / n_tau
    acc = mp.mpf(0)
    for j in range(n_tau + 1):
        tau = -tau_max + j * h
        s = sigma + 1j * tau
        w = mp.mpf("0.5") if j in (0, n_tau) else mp.mpf(1)
        acc += w * abs(F_bump(s, z0)) ** 2
    acc *= h / (2 * PI)
    return mp.sqrt(acc)


def sum_a(lam, kmax: int = 40):
    return sum((a_k(k, lam) for k in range(kmax)), mp.mpf(0))


def psi_trunc_L2(lam, T: float, n_tau: int = 800, kmax: int = 24):
    """Truncated L2 norm of formal Hardy dual psi on [-T,T]."""
    aks = [a_k(k, lam) for k in range(kmax)]
    sks = [mp.mpf("0.5") + k for k in range(kmax)]
    h = mp.mpf(2 * T) / n_tau
    acc = mp.mpf(0)
    for j in range(n_tau + 1):
        tau = -T + j * h
        w = mp.mpf("0.5") if j in (0, n_tau) else mp.mpf(1)
        s_line = SIGMA + 1j * tau
        psi = mp.mpf(0)
        for ak, sk in zip(aks, sks):
            psi += ak / (s_line - sk)
        psi /= 2 * PI
        acc += w * abs(psi) ** 2
    acc *= h
    return mp.sqrt(acc)


def main():
    print("=" * 72)
    print("W1-D1 dual probe (NUMERIC-FORECAST; not a proof gate)")
    print(f"  sigma={SIGMA}  dps={DPS}")
    print("=" * 72)

    # --- coefficient structure ---
    lam0 = mp.mpf("0.3")  # delta = 0.3 + 0.8 = 1.1 at sigma=0.4
    print("\n[P0] even-ladder |a_k| (lam=0.3)")
    for k in range(8):
        ak = a_k(k, lam0)
        print(f"  k={k}  |a|={mp.nstr(abs(ak), 6)}  sign={'+' if ak > 0 else '-'}")
    S_a = sum_a(lam0)
    print(f"  sum_k a_k = {mp.nstr(S_a, 8)}   (psi ~ this/(2 pi i tau) at infinity)")
    print(f"  => C_inf != 0, but 1/tau IS in L2(|tau|>1); check P2 for ||psi|| saturation")

    # --- P1: empirical ratios on bump family ---
    print("\n[P1] R(z0) = |D1(F_z0)| / ||F_z0||_L2(line)   lam=0.3  delta=1.1")
    z0s = [mp.mpf(z) for z in (2, 0, -2, -4, -6, -8)]
    rows = []
    for z0 in z0s:
        d1 = D1_even(lam0, z0)
        nrm = F_L2_line(lam0, z0)
        R = abs(d1) / nrm if nrm > 0 else mp.mpf("nan")
        rows.append((z0, d1, nrm, R))
        print(
            f"  z0={mp.nstr(z0, 4):>6}  |D1|={mp.nstr(abs(d1), 6)}  "
            f"||F||={mp.nstr(nrm, 6)}  R={mp.nstr(R, 6)}"
        )
    R_bulk = rows[0][3]
    R_edge = rows[-1][3]
    print(f"  bulk/edge R ratio (z0=2 vs -8): {mp.nstr(R_bulk / R_edge, 5)}")
    print("  NOTE: edge can make R LARGE (D1 scalar) while full resolvent gain annihilates")
    print("  (OBS-07) — two-stage cancel: G~X[C0-C0(1/2)] small, C0(1/2) need not be.")

    # --- P3: delta scaling at bulk z0=0 ---
    print("\n[P3] delta scaling of R at z0=0  (lam = -2*sigma + delta)")
    for delta in [mp.mpf(d) for d in (1.5, 1.1, 0.6, 0.3, 0.15)]:
        lam = -2 * SIGMA + delta
        d1 = D1_even(lam, mp.mpf(0))
        nrm = F_L2_line(lam, mp.mpf(0))
        R = abs(d1) / nrm
        print(
            f"  delta={mp.nstr(delta, 4)}  lam={mp.nstr(lam, 4)}  "
            f"R={mp.nstr(R, 6)}  R*delta={mp.nstr(R * delta, 6)}"
        )

    # --- P2: truncated dual L2 ---
    print("\n[P2] truncated ||psi||_L2([-T,T]) for formal H^2 dual (lam=0.3)")
    prev = None
    for T in (5.0, 10.0, 20.0, 40.0, 80.0):
        nrm = psi_trunc_L2(lam0, T)
        extra = ""
        if prev is not None:
            # pure 1/tau log growth: || ||^2 ~ log T => || || ~ sqrt(log)
            extra = f"  ratio_to_prev={mp.nstr(nrm / prev, 5)}"
        print(f"  T={T:5.1f}  ||psi||_trunc={mp.nstr(nrm, 7)}{extra}")
        prev = nrm
    print("  Flat ||psi|| => abstract H^2 dual in L2. Separate check: Cauchy recovery")
    print("  of F(s_k) from the line FAILS on bilateral bumps (see findings 2026-08-07).")

    # --- outline flags ---
    print("\n[OUTLINE FLAGS]")
    print("  A. H^2 dual psi in L2, but bilateral L2_sigma is NOT H^2 — Cauchy recovery")
    print("     fails on the repo dense class. A is NOT a D1 discharge for L2_sigma.")
    print("  B. Coboundary telescope: reformulation only.")
    print("  C. Bent-contour Barnes (desk §8): still the named classical route.")
    print("  D. Slit-plane kernel for D(1): best next dual probe.")
    print("  Do NOT discharge via Schur S_off or OBS-07 gains (adversarial 07-25).")
    print("\nDone. No PASS/FAIL gate — see docs/findings_session_2026_08_07_w1_d1_reopen.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
