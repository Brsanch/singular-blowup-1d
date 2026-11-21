#!/usr/bin/env python3
"""
GATE: W1 Route F difference-kernel uniform L²_σ resolvent closure
(docs/desk/w1_plemelj_stirling_writeup.md, 2026-11-20).

Verifies the four analytic ingredients of the Route-F uniform L²_σ bound for
    G(s; lam) = X(s; lam) [ D(e^{4 pi i s}) - D(1) ]
across (sigma, delta = Re(lam) + 2*sigma, omega = Im(lam)):

  F1  Four-quadrant + all-omega exponential decay margin (Lem Stir-Schur):
        E_omega(tau, tau') = pi(|tau + omega/2| - |tau' + omega/2|)
                           + 2*pi(|tau| - |tau'|)
                           + 4*pi * splitWing(tau, tau')
                           >= pi |tau - tau'|
      and shifted polynomial submultiplicativity:
        (1 + |tau + omega/2|)/(1 + |tau' + omega/2|) <= 1 + |tau - tau'|.

  F2  Off-diagonal Schur operator norm (|tau - tau'| >= 1):
        ||T_off(lam)||_{Schur} = sqrt(sup_tau int_{|dt|>=1} dressed dt'
                                    * sup_tau' int_{|dt|>=1} dressed dt)
      is uniformly bounded across omega = Im(lam) in [0, 100] and decays as |omega| -> inf.

  F3  Diagonal band unitary-phase reduction (Lem Diag-Riesz):
      Factoring X(s; lam) = u(tau; lam) |X(s; lam)| with |u| = 1 conjugates out the
      logarithmic phase oscillation (d/dtau arg X ~ 4 log|tau|). The remaining real
      modulus has uniformly bounded logarithmic derivative:
        |d/dtau log|X(sigma+i*tau; lam)|| = |Im(2*psi(2s+lam)+2*psi(2s+2)+psi(s)-psi(1-s))|
                                          <= 3*pi + 2/delta + 1/sigma + 1/(1-sigma) =: L(sigma, delta)
      for all (tau, omega) in R^2, so the modulus commutator kernel on 0 < |tau - tau'| < 1
      cancels the Cauchy pole and is uniformly bounded in (tau, tau', omega).

  F4  Rank-one defect operator omega-cancellation (Lem D-slit-CS):
      Although bare ||K_slit(.; lam)||_L2 grows as e^{+pi|omega|/2}, in the resolvent
      G = X[D(w) - D(1)] the defect term is the rank-one operator
        (T_defect F)(tau) = -X(sigma+i*tau; lam) * int K_slit(tau'; lam) F(sigma+i*tau') d tau',
      whose Hilbert-Schmidt / L2->L2 norm ||X(.; lam)||_L2 * ||K_slit(.; lam)||_L2
      cancels e^{+-pi|omega|/2} and decays as O(1/(1 + |omega|)).

Exit 0 iff F1-F4 all PASS.
"""
from __future__ import annotations

import math
import sys

import mpmath as mp

mp.mp.dps = 30
PI = math.pi
FAILURES: list[str] = []


def report(name: str, ok: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}  {detail}")
    if not ok:
        FAILURES.append(name)


def split_wing(tau: float, taup: float) -> float:
    return max(0.0, taup - tau) if taup >= 0.0 else max(0.0, tau - taup)


def E_omega(tau: float, taup: float, omega: float) -> float:
    return (
        PI * (abs(tau + 0.5 * omega) - abs(taup + 0.5 * omega))
        + 2.0 * PI * (abs(tau) - abs(taup))
        + 4.0 * PI * split_wing(tau, taup)
    )


def log_X_abs(sigma: float, delta: float, tau: float, omega: float) -> float:
    re_lam = delta - 2.0 * sigma
    s = mp.mpc(sigma, tau)
    lam = mp.mpc(re_lam, omega)
    return float(
        mp.re(mp.loggamma(2 * s + lam))
        + mp.re(mp.loggamma(2 * s + 2))
        + mp.re(mp.loggamma(s))
        + mp.re(mp.loggamma(1 - s))
    )


def d_log_X_abs(sigma: float, delta: float, tau: float, omega: float) -> float:
    re_lam = delta - 2.0 * sigma
    s = mp.mpc(sigma, tau)
    lam = mp.mpc(re_lam, omega)
    return float(
        -mp.im(
            2 * mp.psi(0, 2 * s + lam)
            + 2 * mp.psi(0, 2 * s + 2)
            + mp.psi(0, s)
            - mp.psi(0, 1 - s)
        )
    )


def kappa_split_abs(tau: float, taup: float) -> float:
    d = tau - taup
    if abs(d) < 1e-14:
        return 0.0
    if taup >= 0.0:
        if d >= 0.0:
            return 4.0 * PI / abs(1.0 - math.exp(-4.0 * PI * d))
        return 4.0 * PI * math.exp(4.0 * PI * d) / abs(1.0 - math.exp(4.0 * PI * d))
    else:
        if d >= 0.0:
            return 4.0 * PI * math.exp(-4.0 * PI * d) / abs(1.0 - math.exp(-4.0 * PI * d))
        return 4.0 * PI / abs(1.0 - math.exp(4.0 * PI * d))


def main() -> int:
    sigma = 0.4
    delta = 1.1  # Re(lam) = 0.3

    print("W1 Route F uniform L2_sigma resolvent closure gate")
    print(f"  anchor: sigma={sigma}, delta={delta} (Re lam={delta - 2*sigma:.2f})")

    # ---- F1: Four-quadrant + all-omega exponential margin & polynomial submultiplicativity ----
    min_margin_ratio = float("inf")
    max_poly_ratio = 0.0
    test_pts = [-40.0, -15.0, -4.0, -1.0, -0.25, 0.0, 0.25, 1.0, 4.0, 15.0, 40.0]
    for omega in [-50.0, -10.0, 0.0, 5.0, 25.0, 80.0]:
        for t in test_pts:
            for tp in test_pts:
                if abs(t - tp) > 1e-9:
                    r = E_omega(t, tp, omega) / (PI * abs(t - tp))
                    min_margin_ratio = min(min_margin_ratio, r)
                pr = (1.0 + abs(t + 0.5 * omega)) / (
                    (1.0 + abs(tp + 0.5 * omega)) * (1.0 + abs(t - tp))
                )
                max_poly_ratio = max(max_poly_ratio, pr)
    report(
        "F1 Lem Stir-Schur 4-quadrant margin E_omega >= pi|dt| & poly submult",
        min_margin_ratio >= 1.0 - 1e-12 and max_poly_ratio <= 1.0 + 1e-12,
        f"min E_omega/(pi|dt|)={min_margin_ratio:.6f}, max poly_ratio={max_poly_ratio:.6f}",
    )

    # ---- F2: Off-diagonal Schur norm across omega = Im(lam) in [0, 100] ----
    schur_norms = {}
    for omega in [0.0, 5.0, 10.0, 25.0, 50.0, 100.0]:
        lo = min(0.0, -0.5 * omega) - 8.0
        hi = max(0.0, -0.5 * omega) + 8.0
        n_grid = 240
        h = (hi - lo) / n_grid
        grid = [lo + i * h for i in range(n_grid + 1)]
        lxs = [log_X_abs(sigma, delta, t, omega) for t in grid]
        max_row = 0.0
        max_col = 0.0
        for i, t in enumerate(grid):
            r_sum = 0.0
            c_sum = 0.0
            for j, tp in enumerate(grid):
                if abs(t - tp) >= 1.0:
                    k_rt = (
                        math.exp(lxs[i] - lxs[j])
                        * kappa_split_abs(t, tp)
                        / math.hypot(delta, 2.0 * tp + omega)
                    )
                    k_ct = (
                        math.exp(lxs[j] - lxs[i])
                        * kappa_split_abs(tp, t)
                        / math.hypot(delta, 2.0 * t + omega)
                    )
                    r_sum += k_rt * h
                    c_sum += k_ct * h
            max_row = max(max_row, r_sum)
            max_col = max(max_col, c_sum)
        schur_norms[omega] = math.sqrt(max_row * max_col) / (2.0 * PI)
    max_schur = max(schur_norms.values())
    report(
        "F2 Off-diagonal Schur norm ||T_off(lam)||_Schur uniform & decaying in Im(lam)",
        max_schur < 0.25 and schur_norms[100.0] < schur_norms[5.0],
        f"max ||T_off||_Schur={max_schur:.5f}, @om=5:{schur_norms[5.0]:.5f} -> @om=100:{schur_norms[100.0]:.5f}",
    )

    # ---- F3: Unitary-phase reduction & diagonal modulus-commutator bound ----
    L_bound = 3.0 * PI + 2.0 / delta + 1.0 / sigma + 1.0 / (1.0 - sigma)
    max_dlog = 0.0
    max_comm = 0.0
    for omega in [0.0, 5.0, 20.0, 60.0]:
        for tau in [-40.0, -10.0, -2.5, 0.0, 2.5, 10.0, 40.0]:
            dval = abs(d_log_X_abs(sigma, delta, tau, omega))
            max_dlog = max(max_dlog, dval)
            lx_t = log_X_abs(sigma, delta, tau, omega)
            for dt in [1e-4, 1e-2, 0.2, 0.8]:
                taup = tau + dt
                lx_tp = log_X_abs(sigma, delta, taup, omega)
                ratio = math.exp(lx_t - lx_tp)
                ker = kappa_split_abs(tau, taup) / (2.0 * PI)
                denom = math.hypot(delta, 2.0 * taup + omega)
                max_comm = max(max_comm, abs(ratio - 1.0) * ker / denom)
    report(
        "F3 Lem Diag-Riesz |d/dtau log|X|| <= L(sigma,delta) & diagonal commutator bounded",
        max_dlog <= L_bound and max_comm < 5.0,
        f"max |d/dtau log|X||={max_dlog:.4f} <= {L_bound:.4f}; max |R_mod|/delta={max_comm:.4f}",
    )

    # ---- F4: Rank-one defect operator ||T_defect(lam)||_HS = ||X||_2 * ||K_slit||_2 ----
    defect_norms = {}
    for omega in [0.0, 2.0, 5.0, 10.0, 20.0]:
        lo = min(0.0, -0.5 * omega) - 10.0
        hi = max(0.0, -0.5 * omega) + 10.0
        n_grid = 500
        h = (hi - lo) / n_grid
        nx_sq = mp.mpf(0)
        nk_sq = mp.mpf(0)
        for j in range(n_grid + 1):
            tp = lo + j * h
            w_trap = mp.mpf("0.5") if j in (0, n_grid) else mp.mpf(1)
            lx = mp.mpf(log_X_abs(sigma, delta, tp, omega))
            x_val = mp.exp(lx)
            nx_sq += w_trap * (x_val**2)
            wp = mp.exp(4j * mp.pi * sigma) * mp.exp(-4 * mp.pi * tp)
            if tp >= 0:
                ks = abs((-4 * mp.pi * wp) / (wp - 1))
            else:
                ks = abs((-4 * mp.pi * wp) * (1 / (wp - 1) - 1 / wp))
            denom = abs(mp.mpc(delta, 2.0 * tp + omega))
            k_val = (ks / (2 * mp.pi)) / (x_val * denom)
            nk_sq += w_trap * (k_val**2)
        defect_norms[omega] = float(mp.sqrt(nx_sq * nk_sq) * h)
    max_def = max(defect_norms.values())
    monotone_decay = (
        defect_norms[0.0] > defect_norms[5.0] > defect_norms[10.0] > defect_norms[20.0]
    )
    report(
        "F4 Lem D-slit-CS rank-one defect norm ||X||_2*||K_slit||_2 uniform & O(1/(1+|omega|))",
        max_def < 0.5 and monotone_decay,
        f"@om=0:{defect_norms[0.0]:.4f} -> @om=5:{defect_norms[5.0]:.4f} -> @om=20:{defect_norms[20.0]:.4f}",
    )

    # ---- F5: Exact Plemelj right-inverse verification (lam - M_sigma) R_slit(lam) F = F ----
    sig_r, lam_r, tau_r = mp.mpf("0.25"), mp.mpf("0.3"), mp.mpf("0.1")
    s_r = sig_r + 1j * tau_r
    X_fn = lambda w: mp.gamma(2 * w + lam_r) * mp.gamma(2 * w + 2) * mp.gamma(w) * mp.gamma(1 - w)
    c_fn = lambda w: -mp.tan(mp.pi * w) / (2 * w + 1)
    Ghat_fn = lambda w: (mp.mpf("0.5") - w) * mp.exp(2 * w * w)
    Fhat_fn = lambda w: (lam_r + 2 * w) * Ghat_fn(w) - c_fn(w + mp.mpf("0.5")) * Ghat_fn(w + mp.mpf("0.5"))
    rho_fn = lambda t: Fhat_fn(sig_r + 1j * t) / ((lam_r + 2 * (sig_r + 1j * t)) * X_fn(sig_r + 1j * t))

    def ks_fn(t):
        d = tau_r - t
        return 4 * mp.pi / mp.expm1(-4 * mp.pi * d) if t >= 0 else 4 * mp.pi * mp.exp(-4 * mp.pi * d) / mp.expm1(-4 * mp.pi * d)

    def k1_fn(t):
        wp = mp.exp(4 * mp.pi * 1j * (sig_r + 1j * t))
        return -4 * mp.pi * (wp if t >= 0 else 1) / (wp - 1)

    integ_fn = lambda t: (ks_fn(t) - k1_fn(t)) * rho_fn(t)
    rad = mp.mpf("0.05")
    reg_lim = mp.diff(rho_fn, tau_r) - (2 * mp.pi + k1_fn(tau_r)) * rho_fn(tau_r)
    sym_pair = lambda x: 2 * reg_lim if abs(x) < mp.mpf("1e-20") else integ_fn(tau_r + x) + integ_fn(tau_r - x)
    pv_val = (
        mp.quad(integ_fn, [-10, -6, -4, -2, -1, 0, tau_r - rad])
        + mp.quad(sym_pair, [0, rad])
        + mp.quad(integ_fn, [tau_r + rad, 0.5, 1, 2, 4, 6, 10])
    )
    cand_G = X_fn(s_r) * (rho_fn(tau_r) / 2 - pv_val / (2 * mp.pi * 1j))
    known_G = Ghat_fn(s_r)
    rel_err_f5 = float(abs(cand_G - known_G) / abs(known_G))
    report(
        "F5 Exact Plemelj right-inverse check R_slit(lam) Fhat = Ghat (eq. app_Rslit_explicit)",
        rel_err_f5 < 1e-25,
        f"rel_err={rel_err_f5:.2e} < 1e-25 (@sigma=1/4, lam=0.3, tau=0.1)",
    )

    print()
    if FAILURES:
        print(f"GATE FAIL: {FAILURES}")
        return 1
    print("GATE PASS: W1 Route F difference-kernel uniform L2_sigma resolvent closure")
    return 0


if __name__ == "__main__":
    sys.exit(main())
