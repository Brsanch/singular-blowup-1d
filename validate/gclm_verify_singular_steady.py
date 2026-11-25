#!/usr/bin/env python3
"""
GATE 1 (Target B -- gCLM singular profiles, arXiv:2603.25104v4):
Independently certify the exact singular self-similar steady states for a <= 0.

Checks:
  1. For a grid of a < 0 spanning all three regimes (-1 < a < 0, a = -1, a < -1):
     a. Independent adaptive-quadrature Hilbert transform H(Omega_bar_a)(X) matches
        analytic Ux_bar(X) on both X < 1 and X > 1 to < 1e-11.
     b. Analytic U_bar_a(0) = 0, H(Omega_bar_a)(0) = 1, and exact edge stagnation
        c_l_bar * 1 + a * U_bar_a(1) = 0 to machine precision.
     c. Steady residual of dynamic-rescaling equation (2.7) on X > 1 is machine zero
        (< 1e-13).
     d. Background advection formula A_bar_a(1+v) = (1-a)*(v - cos(pi*mu)*v^{1-mu})
        and regime classification:
          -1 < a < 0: strictly positive on X > 1 (finite-time inflow at X=1+);
          a = -1:     linear 2*(X-1) (logarithmic inflow, Hou-Luo transport);
          a < -1:     interior repelling stagnation point X_* = 1 + (cos(pi*mu))^{1/mu} > 1.
  2. Exact CLM (a = 0) singular limits (Theorem 1 odd & Corollary 2 half-line):
     a. Direct dynamic evolution of h = F + i G from degenerate smooth initial data
        converges pointwise away from X=1 (and X=+-1) to H(Omega_bar) = 2/(1-X) (half-line)
        and 2/(1-X^2) (odd), and weakly to -2*pi*delta(X-1) and -pi*(delta(X-1)-delta(X+1)).
  3. Exact Mellin symbol c_a(s) verification on test power functions G(1+v) = v^alpha,
     plus exact amplitude zero mode c_a(mu) = 0 and exact a = -1 conjugation identity
     N_{-1}(s) * c_{-1}(s) = -N_{-1}(s - 1/2).
"""

import os
import sys
import numpy as np
from scipy.integrate import quad

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from gclm_singular.profiles import (  # noqa: E402
    mu_of_a,
    scaling_exponents,
    omega_bar,
    Ux_bar,
    U_bar,
    advection_bar,
    steady_residual,
    interior_stagnation_point,
    hilbert_omega_bar_quad,
    mellin_symbol_c,
    conjugation_N_am1,
)


def check_a_negative():
    a_grid = [-2.5, -2.0, -1.5, -1.0, -0.75, -0.5, -0.25, -0.1]
    X_quad = np.array([-3.0, -0.5, 0.0, 0.5, 0.9, 1.1, 1.5, 2.5, 10.0, 50.0])
    X_res = 1.0 + np.geomspace(1e-4, 1e3, 2000)

    max_quad_err = 0.0
    max_stag_err = 0.0
    max_res_err = 0.0
    max_adv_err = 0.0
    regime_ok = True

    print("  [1] Explicit singular steady states (Theorem 3, a < 0):")
    print(f"      {'a':>6s} {'mu':>8s} {'c_l':>6s} {'X_*':>10s} "
          f"{'max|H_q-Ux|':>12s} {'|A(1)|':>10s} {'max|R_ss|':>11s}")

    for a in a_grid:
        mu = mu_of_a(a)
        c_l, c_w, gamma = scaling_exponents(a)

        # (a) Quadrature Hilbert vs analytic Ux_bar
        Hq = hilbert_omega_bar_quad(X_quad, a)
        Ux_an = Ux_bar(X_quad, a)
        err_q = float(np.max(np.abs(Hq - Ux_an)))
        max_quad_err = max(max_quad_err, err_q)

        # (b) Boundary/stagnation checks
        u0_err = abs(float(U_bar(0.0, a)))
        ux0_err = abs(float(Ux_bar(0.0, a)) - 1.0)
        stag_err = abs(float(advection_bar(1.0, a)))
        max_stag_err = max(max_stag_err, u0_err, ux0_err, stag_err)

        # (c) Steady residual relative to term-sum scale (no cancellation-masking)
        R = steady_residual(X_res, a)
        Om = omega_bar(X_res, a)
        OmX = mu * np.sin(np.pi * mu) * (X_res - 1.0) ** (-1.0 - mu)
        Ux = Ux_bar(X_res, a)
        U = U_bar(X_res, a)
        term_scale = (
            np.abs(c_l * X_res * OmX)
            + np.abs(a * U * OmX)
            + np.abs(c_w * Om)
            + np.abs(Ux * Om)
        )
        rel_R = float(np.max(np.abs(R) / term_scale))
        max_res_err = max(max_res_err, rel_R)

        # (d) Advection identity A_bar(1+v) = (1-a)*(v - cos(pi*mu)*v^{1-mu})
        v = X_res - 1.0
        A_direct = advection_bar(X_res, a)
        A_formula = (1.0 - a) * (v - np.cos(np.pi * mu) * (v ** (1.0 - mu)))
        adv_err = float(np.max(np.abs(A_direct - A_formula) / (1.0 + np.abs(A_formula))))
        max_adv_err = max(max_adv_err, adv_err)

        # Regime verification
        x_star = interior_stagnation_point(a)
        if a > -1.0:
            regime_ok = regime_ok and (x_star is None) and bool(np.all(A_direct > 0.0))
            x_star_str = "None (>0)"
        elif a == -1.0:
            regime_ok = regime_ok and (x_star == 1.0) and bool(np.all(A_direct > 0.0))
            x_star_str = "1.0 (edge)"
        else:
            A_star = abs(float(advection_bar(x_star, a)))
            A_left = float(advection_bar(1.0 + 0.5 * (x_star - 1.0), a))
            A_right = float(advection_bar(x_star + 0.5, a))
            regime_ok = regime_ok and (A_star < 1e-13) and (A_left < 0.0) and (A_right > 0.0)
            x_star_str = f"{x_star:.5f}"

        print(f"      {a:+6.2f} {mu:8.5f} {c_l:6.2f} {x_star_str:>10s} "
              f"{err_q:12.3e} {stag_err:10.3e} {rel_R:11.3e}")

    print(f"      max advection identity rel err = {max_adv_err:.3e}, regime_ok = {regime_ok}")
    ok = (
        max_quad_err < 1e-11
        and max_stag_err < 1e-14
        and max_res_err < 1e-13
        and max_adv_err < 1e-14
        and regime_ok
    )
    return ok


def check_a_zero_clm():
    """Verify Theorem 1 (odd) and Corollary 2 (half-line) exact CLM singular convergence."""
    # Odd degenerate initial datum: F0(x) = -x^3 / (1 + x^2)^2? Need G0 = H(F0) with G0(0)=2, G0''(0)=4.
    # Exact holomorphic function in UHP: h0(z) = 2i / (1 - i z)^2? Wait:
    # Let h0(z) = i * (2 + 2 z^2) / (1 - i z)^4? Simpler:
    # Use a rational odd F0 with exact Hilbert transform!
    # Recall for z in UHP, f(z) = c / (z + i)^k is holomorphic -> Re = F, Im = H(F).
    # For odd F0, take pure imaginary Laurent at z = i:
    #   h(z) = i * [ A / (1 + z^2) + B / (1 + z^2)^2 ]? No, poles at +-i are NOT UHP-holomorphic.
    # UHP-holomorphic pole ONLY at z = -i:
    #   (z + i)^{-1} = (z - i)/(z^2 + 1) -> F = z/(z^2+1), G = -1/(z^2+1) (non-degenerate, F'(0)=1).
    # Take linear combination cancelling z^1 term in Re h0(x):
    #   Let phi_k(z) = (i / (z + i))^k. On R:
    #     k=1: i/(x+i) = (1 + i x)/(1+x^2) -> Re = 1/(1+x^2) (even), Im = x/(1+x^2) (odd).
    #   Multiply by i so Re is odd and Im is even!
    #     psi_1(z) = i * i/(z+i) = -1/(z+i) -> Re = -x/(1+x^2), Im = 1/(1+x^2).
    #     psi_2(z) = i * (i/(z+i))^2 = -i/(z+i)^2 = -i(x^2-1 - 2ix)/(1+x^2)^2
    #              -> Re = -2x/(1+x^2)^2, Im = (1 - x^2)/(1+x^2)^2.
    #   Wait: to get Re(h0) = O(x^3) and Re(h0) < 0 for x > 0:
    #     psi_1: Re = -x/(1+x^2) = -x + x^3 - ...
    #     psi_2: Re = -2x/(1+x^2)^2 = -2x + 4x^3 - ...
    #     2*psi_1 - psi_2: Re = -2x/(1+x^2) + 2x/(1+x^2)^2 = -2 x^3 / (1 + x^2)^2  (< 0 for x > 0!!),
    #                      Im = 2/(1+x^2) - (1-x^2)/(1+x^2)^2 = (1 + 3 x^2) / (1 + x^2)^2!
    #   Check Im at x=0: Im(0) = 1, Im''(0) = 2*(3 - 2) = 2!
    #   So scale by alpha = 2: h0(x) = 2 * (2*psi_1(x) - psi_2(x)):
    #     F0(x) = -4 x^3 / (1 + x^2)^2  (odd, <0 on x>0, F0'(0)=F0''(0)=0!),
    #     G0(x) = 2*(1 + 3 x^2) / (1 + x^2)^2  -> G0(0) = 2, G0''(0) = 4! EXACTLY!
    def h0_odd(x):
        F0 = -4.0 * (x ** 3) / ((1.0 + x ** 2) ** 2)
        G0 = 2.0 * (1.0 + 3.0 * x ** 2) / ((1.0 + x ** 2) ** 2)
        return F0 + 1j * G0

    # Exact dynamic rescaling solution (Section 3, eq. line 604) at tau:
    tau = 20.0
    xs = np.array([0.2, 0.5, 0.8, 1.3, 2.0, 5.0])
    h_tau = h0_odd(xs * np.exp(-0.5 * tau)) / (
        np.exp(tau) + 0.5j * (np.exp(tau) - 1.0) * h0_odd(xs * np.exp(-0.5 * tau))
    )
    F_err_odd = float(np.max(np.abs(h_tau.real)))
    G_exact_odd = 2.0 / (1.0 - xs ** 2)
    G_err_odd = float(np.max(np.abs(h_tau.imag - G_exact_odd)))

    # Cancellation-free evaluation of F(x, tau) using 2 - G0(y) = 2*(y^4 - y^2)/(1+y^2)^2
    # (Manuscript.tex line 638). Allows clean quadrature at tau_w = 18.
    tau_w = 18.0
    et = np.exp(tau_w)
    em = et - 1.0
    eh = np.exp(-0.5 * tau_w)

    def F_odd_at(x):
        y = x * eh
        y2 = y * y
        d2 = (1.0 + y2) ** 2
        F0 = -4.0 * (y * y2) / d2
        two_minus_G0 = 2.0 * (y2 * y2 - y2) / d2
        denom_re = 2.0 + em * two_minus_G0
        return float(4.0 * et * F0 / ((em * F0) ** 2 + denom_re ** 2))

    phi = lambda x: np.exp(-((x - 1.0) ** 2))
    # Peak is at x = 1 + O(e^{-tau_w/2}) with width O(e^{-tau_w/2}) ~ 1.2e-4
    w = 1.0 / np.sqrt(et)
    pts = sorted({0.0, 0.5, 0.99, 1.0 - 50 * w, 1.0 - 5 * w, 1.0, 1.0 + 5 * w, 1.0 + 50 * w, 1.01, 2.5})
    weak_odd = 0.0
    for lo, hi in zip(pts[:-1], pts[1:]):
        val, _ = quad(lambda x: F_odd_at(x) * phi(x), lo, hi, epsabs=1e-11, epsrel=1e-10, limit=200)
        weak_odd += val
    weak_odd_err = abs(weak_odd - (-np.pi * phi(1.0)))

    # Half-line exact steady check (Corollary 2): via X -> sqrt(X) transformation from odd case
    Xs_hl = xs ** 2
    G_exact_hl = 2.0 / (1.0 - Xs_hl)
    G_err_hl = float(np.max(np.abs(h_tau.imag - G_exact_hl)))

    weak_hl = 0.0
    for lo, hi in zip(pts[:-1], pts[1:]):
        val, _ = quad(
            lambda y: F_odd_at(y) * phi(y ** 2) * (2.0 * y),
            lo, hi, epsabs=1e-11, epsrel=1e-10, limit=200,
        )
        weak_hl += val
    weak_hl_err = abs(weak_hl - (-2.0 * np.pi * phi(1.0)))

    print("\n  [2] Exact CLM (a = 0) singular limits (Theorem 1 odd & Corollary 2 half-line):")
    print(f"      odd  (tau={tau:.0f}): max|F|={F_err_odd:.3e}, max|G - 2/(1-X^2)|={G_err_odd:.3e}")
    print(f"      odd  weak (tau={tau_w:.0f}): |int F*phi - (-pi*phi(1))| = {weak_odd_err:.3e}")
    print(f"      half (tau={tau:.0f}): max|G - 2/(1-X)|  ={G_err_hl:.3e}")
    print(f"      half weak (tau={tau_w:.0f}): |int F*phi - (-2pi*phi(1))|= {weak_hl_err:.3e}")

    return (
        F_err_odd < 1e-3
        and G_err_odd < 1e-6
        and weak_odd_err < 1e-3
        and G_err_hl < 1e-6
        and weak_hl_err < 2e-3
    )


def check_mellin_symbol_and_conjugation():
    """Verify exact Mellin symbol c_a(s) on power functions v^alpha and a=-1 conjugation."""
    a_grid = [-2.0, -1.5, -1.0, -0.5, -0.2]
    max_pow_err = 0.0
    max_zero_err = 0.0

    for a in a_grid:
        mu = mu_of_a(a)
        # Amplitude mode s = mu -> c_a(mu) = 0
        max_zero_err = max(max_zero_err, abs(mellin_symbol_c(mu, a)))

        # Test G(1+v) = v^alpha with 0 < mu - alpha < 1 (so s = mu - alpha in (0, 1))
        for s_test in [0.15, 0.35, 0.65, 0.85]:
            alpha = mu - s_test
            # Direct physical-space action on G(1+v) = v^alpha:
            #   g(1+v) = -sin(pi*mu) * v^{-s_test}
            #   H_+[g](1+v) = +sin(pi*mu)*cot(pi*s_test)*v^{-s_test}
            #   Q_+[g](1+v) = sin(pi*mu)*cot(pi*s_test)/(1-s_test) * v^{-s_test}
            #   -A_bar * G_v = -(1-a)*alpha * v^alpha + (1-a)*cos(pi*mu)*alpha * v^{alpha-mu}
            # Total coefficient of v^{alpha - mu} (= v^{-s_test}):
            coeff_phys = (
                (1.0 - a) * np.cos(np.pi * mu) * alpha
                + np.sin(np.pi * mu) * (np.cos(np.pi * s_test) / np.sin(np.pi * s_test))
                * (1.0 + a * mu / (1.0 - s_test))
            )
            coeff_mellin = float(np.real(mellin_symbol_c(s_test, a)))
            max_pow_err = max(max_pow_err, abs(coeff_phys - coeff_mellin))

    # Check exact a = -1 conjugation: N_{-1}(s) * c_{-1}(s) + N_{-1}(s - 1/2) == 0 (relative)
    s_complex = np.array([0.2 + 0.3j, 0.35 - 1.7j, 0.1 + 4.2j, 0.45 + 0.0j])
    rhs_conj = conjugation_N_am1(s_complex - 0.5)
    conj_res = conjugation_N_am1(s_complex) * mellin_symbol_c(s_complex, -1.0) + rhs_conj
    max_conj_err = float(np.max(np.abs(conj_res) / np.abs(rhs_conj)))

    print("\n  [3] Exact Mellin symbol c_a(s) & a=-1 conjugation:")
    print(f"      max |c_a(mu)| (amplitude zero mode)           = {max_zero_err:.3e}")
    print(f"      max |coeff_phys - c_a(mu-alpha)| (power test) = {max_pow_err:.3e}")
    print(f"      max rel |N_-1(s) c_-1(s) + N_-1(s-1/2)|       = {max_conj_err:.3e}")

    return max_zero_err < 1e-14 and max_pow_err < 1e-13 and max_conj_err < 1e-14


def main():
    print("=" * 78)
    print("  GATE 1: gCLM SINGULAR SELF-SIMILAR STEADY STATES (arXiv:2603.25104v4)")
    print("=" * 78)
    ok1 = check_a_negative()
    ok2 = check_a_zero_clm()
    ok3 = check_mellin_symbol_and_conjugation()
    verdict = ok1 and ok2 and ok3
    print("\n" + "-" * 78)
    print(f"  [1] a < 0 steady & quadrature: {ok1} | "
          f"[2] a = 0 CLM limits: {ok2} | "
          f"[3] Mellin & conjugation: {ok3}")
    print("\n  VERDICT:", "PASS -- Gate 1 certified" if verdict else "FAIL")
    print("=" * 78)
    return 0 if verdict else 1


if __name__ == "__main__":
    sys.exit(main())
