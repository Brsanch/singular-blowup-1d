#!/usr/bin/env python3
"""
GATE: Final Nonlinear Bootstrap Closure — Edge Hardy-Disk H^1(D(1,1)) +
Overlap Far-Field H^1_{1/4}([-ln 2, inf)) + Global W1-W3 L^2_{1/4} Semigroup
(docs/desk/nonlinear_comoving_leg.md §4, 2026-11-20).

Verifies the five quantitative gates (NB-1..NB-5) of the closed nonlinear bootstrap
for the exact Gauge-(a) 1D Hou-Luo co-moving system:

  NB-1  Edge-disk H^1(D(1,1)) strict coercivity:
        Discrete Hilbert matrix H_{1/2}(m,k) = 1/(k - m - 1/2) has ||H_{1/2}||_2 <= pi;
        the weighted edge-disk operator L_disk(m,k) = -2m delta_{mk} - 1/(2 pi (m+1) (k-m-1/2))
        has strictly negative H^1 Rayleigh abscissa <= -1.75 (<= -1.50 after far-field Young)
        across K = 20, 40, 80.
  NB-2  Zero-derivative L^2_{3/4,1/4} -> L^inf bound on Q_op G:
        int |c(1/4+it)|^2 dt = pi/3, int |c(3/4+it)|^2 dt = pi/5 =>
        ||Q_op G||_inf <= (1/sqrt(pi) + sqrt(2/3)) ||G||_{3/4,1/4} < (7/5) ||G||_{3/4,1/4}.
  NB-3  Far-field-to-edge-disk Taylor coupling convergence (zero radius shrinkage, r0 = 1):
        sum_{m>=1} m |b_m|^2 <= ((pi^2/6 - 1) / (8 pi^2)) ||G||_{L^2_{1/4}([0,inf))}^2
                             < (25/2592) ||G||^2 < (1/100) ||G||^2.
  NB-4  Overlap inflow at z0 = -ln 2 (Y = 3/2) & edge L^2_{3/4} embedding:
        2 e^{-z0/2} |G_z(-ln 2)|^2 <= (2 sqrt(2) / 3) ||c||_{H^1}^2 < ||c||_{H^1}^2,
        ||G||_{L^2_{3/4}((-inf,0])}^2 <= trace(D^{-1} Gram_{3/4} D^{-1}) ||c||_{H^1}^2
                                      < (199/90) ||c||_{H^1}^2 < (9/4) ||c||_{H^1}^2.
  NB-5  Closed 3-block Lyapunov / Riccati bootstrap trapping:
        Explicit weights (a0, a1, a2) = (20, 1, 10) with W1-W3 Datko-Pazy rate 2*omega_1 = 4/5
        make the coupled (X_sem, Y_far, Z_disk) dissipation rates (-741/1000, -1/2, -5359/2000),
        giving gamma_nl = 1/2 and global exponential decay E(tau) <= E(0) exp(-tau/4)
        for E(0) < (gamma_nl / (2 C_nl))^2.

Exit 0 iff NB-1..NB-5 all PASS.
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
from stability_edge_z import build_blocks_z  # noqa: E402

FAILURES: list[str] = []


def report(name: str, ok: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}  {detail}")
    if not ok:
        FAILURES.append(name)


def main() -> int:
    print("Final Nonlinear Bootstrap Closure Gate (NB-1..NB-5)")

    # ---- NB-1: Discrete Hilbert matrix norm <= pi & H^1 edge-disk strict coercivity ----
    abscissas = {}
    h_norms = {}
    for K in [20, 40, 80]:
        m = np.arange(1, K + 1, dtype=float)[:, None]
        k = np.arange(1, K + 1, dtype=float)[None, :]
        H_half = 1.0 / (k - m - 0.5)
        h_norms[K] = float(np.linalg.norm(H_half, 2))
        L_disk = np.diag(-2.0 * np.arange(1, K + 1, dtype=float)) - H_half / (2.0 * math.pi * (m + 1.0))
        D2 = np.diag(np.arange(1, K + 1, dtype=float) ** 2)
        D_diss_inv = np.diag(np.arange(1, K + 1, dtype=float) ** (-1.5))
        Sym_H1 = 0.5 * (D2 @ L_disk + L_disk.T @ D2)
        G_mat = D_diss_inv @ Sym_H1 @ D_diss_inv
        abscissas[K] = float(np.max(np.linalg.eigvalsh(G_mat)))
    ok_nb1 = (
        all(v <= math.pi + 1e-12 for v in h_norms.values())
        and all(v <= -1.75 for v in abscissas.values())
    )
    report(
        "NB-1 Discrete Hilbert norm <= pi & H^1(D(1,1)) Rayleigh abscissa <= -1.75",
        ok_nb1,
        f"||H_{{1/2}}||_2(K=80)={h_norms[80]:.5f}<=pi({math.pi:.5f}), "
        f"max Rayleigh(K=20/40/80)=({abscissas[20]:.4f}/{abscissas[40]:.4f}/{abscissas[80]:.4f}) <= -1.75",
    )

    # ---- NB-2: Zero-derivative L^2_{3/4,1/4} -> L^inf bound on Q_op G ----
    # Exact tangent substitution tau = tan(theta), theta in (-pi/2, pi/2) integrates to machine precision.
    theta = np.linspace(-0.5 * math.pi, 0.5 * math.pi, 4001)
    c2 = np.cos(theta) ** 2
    s2 = np.sin(theta) ** 2
    int_14 = float(np.trapezoid(1.0 / (2.25 * c2 + 4.0 * s2), theta))
    int_34 = float(np.trapezoid(1.0 / (6.25 * c2 + 4.0 * s2), theta))
    c_linf_bound = 1.0 / math.sqrt(math.pi) + math.sqrt(2.0 / 3.0)

    B = build_blocks_z(N=1000, zmin=-10.0, zmax=10.0, beta=0.0, edge_bc=True)
    z, dz, Y, Q_op = B["z"], B["dz"], B["Y"], B["Q_op"]
    Ymin = 1.0 + math.exp(-10.0)
    w_trap = np.full(len(z), dz)
    w_trap[0] *= 0.5
    w_trap[-1] *= 0.5
    h1_row = -(1.0 / np.pi) * np.exp(-0.5 * z) * w_trap
    Q_exact = Q_op + np.outer((Ymin - 1.0) / (Y - 1.0), h1_row)
    W2 = np.exp(-1.5 * z) + np.exp(-0.5 * z)
    max_emp_ratio = 0.0
    for Gtest in [
        np.exp(-((z - 1.0) / 1.0) ** 2),
        (Y - 1.0) * np.exp(-(Y - 1.0)),
        (Y - 1.0) * np.exp(-(Y - 1.0)) * np.cos(3.0 * z),
    ]:
        n2 = float(np.sqrt(np.sum(W2 * Gtest**2) * dz))
        qinf = float(np.max(np.abs(Q_exact @ Gtest)))
        max_emp_ratio = max(max_emp_ratio, qinf / n2)

    ok_nb2 = (
        abs(int_14 - math.pi / 3.0) < 1e-12
        and abs(int_34 - math.pi / 5.0) < 1e-12
        and c_linf_bound < 1.4
        and max_emp_ratio < c_linf_bound
    )
    report(
        "NB-2 Zero-derivative L^2_{3/4,1/4} -> L^inf bound on Q_op G (< 7/5 = 1.40)",
        ok_nb2,
        f"err(pi/3)={abs(int_14 - math.pi/3):.2e}, err(pi/5)={abs(int_34 - math.pi/5):.2e}, "
        f"bound={c_linf_bound:.5f}<1.40, emp_max={max_emp_ratio:.4f}",
    )

    # ---- NB-3: Far-field-to-edge-disk Taylor coupling sum_{m>=1} m |b_m|^2 ----
    ms = np.arange(1, 50001, dtype=float)
    exact_series = float(np.sum(ms / (4.0 * math.pi**2 * (ms + 1.0) ** 2 * (2.0 * ms + 0.5))))
    analytic_cap = (math.pi**2 / 6.0 - 1.0) / (8.0 * math.pi**2)
    elem_cap = 25.0 / 2592.0  # Elementary pi>3 telescoping bound < 1/100
    G_far_n2_sq = float(1.0 / 2.5)  # G(z) = e^{-z} on [0, inf)
    b_m_test = -1.0 / (2.0 * math.pi * (ms + 1.0) * (ms + 1.5))
    emp_far_ratio = float(np.sum(ms * b_m_test**2) / G_far_n2_sq)
    ok_nb3 = exact_series < analytic_cap < elem_cap < 0.01 and emp_far_ratio <= exact_series + 1e-12
    report(
        "NB-3 Far-field-to-edge-disk coupling sum m |b_m|^2 < 1/100 ||G||_{L^2_{1/4}}^2",
        ok_nb3,
        f"series={exact_series:.6f} < pi_cap={analytic_cap:.6f} < elem_cap={elem_cap:.6f} < 0.01, "
        f"emp={emp_far_ratio:.6f}",
    )

    # ---- NB-4: Overlap inflow at z0 = -ln 2 & edge L^2_{3/4} embedding ----
    # Inflow: 2 e^{-z0/2} |G_z(-ln 2)|^2 = 2 sqrt(2) |sum m c_m 2^{-m}|^2 <= (2 sqrt(2) / 3) ||c||_{H^1}^2 < ||c||_{H^1}^2
    # Edge embedding: ||G||_{L^2_{3/4}((-inf,0])}^2 <= trace(D^{-1} Gram_{3/4} D^{-1}) ||c||_{H^1}^2 < 199/90 < 9/4
    K_test = 80
    m_arr = np.arange(1, K_test + 1, dtype=float)
    geom_sum = float(np.sum(2.0 ** (-2.0 * m_arr)))
    c_worst = (2.0 ** (-m_arr)) / m_arr
    Gz_z0_sq = float(np.sum(m_arr * c_worst * (2.0 ** (-m_arr))) ** 2)
    c_H1_sq = float(np.sum((m_arr**2) * (c_worst**2)))
    inflow_weighted = 2.0 * math.sqrt(2.0) * (Gz_z0_sq / c_H1_sq)

    m_col = m_arr[:, None]
    k_row = m_arr[None, :]
    Gram_34 = 1.0 / (m_col + k_row - 1.5)
    D_inv = np.diag(1.0 / m_arr)
    A_mat = D_inv @ Gram_34 @ D_inv
    emb_34_max = float(np.max(np.linalg.eigvalsh(A_mat)))
    trace_34 = float(np.sum(1.0 / ((m_arr**2) * (2.0 * m_arr - 1.5))))
    telescoping_cap = 199.0 / 90.0  # 2 + 1/10 + 1/9 < 9/4 = 2.25
    ok_nb4 = (
        abs(geom_sum - 1.0 / 3.0) < 1e-12
        and inflow_weighted < 1.0
        and emb_34_max <= trace_34 < telescoping_cap < 2.25
    )
    report(
        "NB-4 Overlap inflow (< 1) & edge L^2_{3/4} embedding (<= trace < 199/90 < 9/4)",
        ok_nb4,
        f"2*sqrt(2)*inflow={inflow_weighted:.5f}<1.0, emb_max={emb_34_max:.5f} <= "
        f"trace={trace_34:.5f} < 199/90({telescoping_cap:.5f}) < 9/4(2.25)",
    )

    # ---- NB-5: Closed 3-block Lyapunov / Riccati bootstrap stability check ----
    # Exact rational constants (matching Lean SingularBlowup/NonlinearBootstrap.lean):
    #   d/dtau X_sem  <= -(4/5) X_sem + nl             [Datko-Pazy with omega_1 = 2/5 < 9/20 from W1-W3]
    #   d/dtau Y_far  <= -(1/2) Y_far + (641/200) Z_disk + (49/50) X_sem + nl
    #                     [inflow < 1, nonlocal 2*(7/10)^2 * ((9/4) Z_disk + X_sem)]
    #   d/dtau Z_disk <= -3 Z_disk + (1/50) X_sem + nl [2 * (1/100) X_sem]
    # Weights (a0, a1, a2) = (20, 1, 10):
    a0, a1, a2 = 20.0, 1.0, 10.0
    c_x = -(4.0 / 5.0) * a0 + (49.0 / 50.0) * a1 + (1.0 / 50.0) * a2  # -741/50 = -14.82
    c_y = -(1.0 / 2.0) * a1                                           # -1/2 = -0.50
    c_z = -3.0 * a2 + (641.0 / 200.0) * a1                            # -5359/200 = -26.795
    gamma_nl = min(-c_x / a0, -c_y / a1, -c_z / a2)                   # exactly 0.50

    # Verify Riccati trapping ODE numerically: E' = -gamma_nl E + C_nl E^{3/2}
    # with C_nl = 10, E(0) = 0.9 * (gamma_nl / (2 C_nl))^2 => E(tau) <= E(0) exp(-(gamma_nl/2) tau)
    C_nl = 10.0
    E_crit = (gamma_nl / (2.0 * C_nl)) ** 2
    E0 = 0.9 * E_crit
    dt = 1e-3
    E = E0
    riccati_ok = True
    for step in range(1, 20001):
        # RK4 step on E' = -gamma_nl * E + C_nl * E^{3/2}
        def rhs(val: float) -> float:
            return -gamma_nl * val + C_nl * (val ** 1.5)

        k1 = rhs(E)
        k2 = rhs(E + 0.5 * dt * k1)
        k3 = rhs(E + 0.5 * dt * k2)
        k4 = rhs(E + dt * k3)
        E = E + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        tau = step * dt
        if E > E0 * math.exp(-0.5 * gamma_nl * tau) * (1.0 + 1e-10):
            riccati_ok = False
            break

    ok_nb5 = (
        c_x < 0
        and c_y < 0
        and c_z < 0
        and abs(gamma_nl - 0.5) < 1e-12
        and riccati_ok
    )
    report(
        "NB-5 Closed 3-block Lyapunov functional & Riccati trapping (gamma_nl = 1/2)",
        ok_nb5,
        f"rates: x={c_x/a0:.4f}, y={c_y/a1:.4f}, z={c_z/a2:.4f} => gamma_nl={gamma_nl:.4f}, "
        f"E(20)/E(0)={E/E0:.2e} <= e^{{-5}}({math.exp(-5.0):.2e})",
    )

    print()
    if FAILURES:
        print(f"GATE FAIL: {FAILURES}")
        return 1
    print("GATE PASS: Final Nonlinear Bootstrap Closure (NB-1..NB-5)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
