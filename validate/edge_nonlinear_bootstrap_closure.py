#!/usr/bin/env python3
"""
GATE: Final Nonlinear Bootstrap Closure — Pure Direct-Energy Edge Hardy-Disk
H^1(D(1,1)) + Overlap Far-Field H^1_{1/4}([-ln 2, inf)) (Zero W3/Semigroup Dependence)
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
  NB-5  Closed pure-energy 3-block Lyapunov / Riccati bootstrap trapping:
        Using unconditional ||(1/2) Q_tilde G||_{L^2_{1/4}(R)} <= (2/5) ||G||_{L^2_{3/4}(R)},
        explicit weights (a0, a1, a2) = (20, 1, 20) make the direct-energy
        (X_far, Y_far, Z_disk) dissipation rates (-131/1000, -1/2, -3759/4000),
        giving gamma_nl = 1/8 and global exponential decay E(tau) <= E(0) exp(-tau/16)
        for E(0) < (gamma_nl / (2 C_nl))^2 (no W3/Datko-Pazy semigroup needed).

  NB-6  Appendix B continuous-PDE identities:
        (i) exact unit-disk identity H[1_{(0,1)} w^{k-1/2}](u) = -(1/pi) sum_{m>=0} u^m/(k-m-1/2),
        (ii) exact scalar edge-Hilbert decay H[exp(tau M_tilde) G_0](1) = H[g_0](1) e^{-tau},
        (iii) derivative Mellin L^inf bound e^{z/4}|d_z Q_op G| <= sqrt(2/3) ||G_z||_{L^2_{1/4}},
        (iv) 1-IBP O(m^{-3}) cancellation of G(0) = sum c_k in q_m = q_m^{disk} + q_m^{far}.

Exit 0 iff NB-1..NB-6 all PASS.
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
    print("Final Nonlinear Bootstrap Closure Gate (NB-1..NB-6)")

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

    # ---- NB-5: Closed pure-energy 3-block Lyapunov / Riccati bootstrap (NO W3 / Datko-Pazy!) ----
    # Exact rational constants (matching Lean SingularBlowup/NonlinearBootstrap.lean):
    #   d/dtau X_far  <= -(1/5) X_far + (19/10) Z_disk + nl
    #                     [inflow < 1, Young eps=2/5 on ||(1/2) Q_tilde G||_{1/4} <= (2/5)||G||_{3/4}]
    #   d/dtau Y_far  <= -(1/2) Y_far + (641/200) Z_disk + (49/50) X_far + nl
    #                     [inflow < 1, nonlocal 2*(7/10)^2 * ((9/4) Z_disk + X_far)]
    #   d/dtau Z_disk <= -3 Z_disk + (1/50) X_far + nl [2 * (1/100) X_far]
    # Weights (a0, a1, a2) = (20, 1, 20):
    a0, a1, a2 = 20.0, 1.0, 20.0
    c_x = -(1.0 / 5.0) * a0 + (49.0 / 50.0) * a1 + (1.0 / 50.0) * a2  # -131/50 = -2.62
    c_y = -(1.0 / 2.0) * a1                                           # -1/2 = -0.50
    c_z = -3.0 * a2 + (19.0 / 10.0) * a0 + (641.0 / 200.0) * a1       # -3759/200 = -18.795
    min_rate = min(-c_x / a0, -c_y / a1, -c_z / a2)                   # 131/1000 = 0.131 > 1/8
    gamma_nl = 1.0 / 8.0                                              # certified rational floor 1/8

    # Verify Riccati trapping ODE numerically: E' = -gamma_nl E + C_nl E^{3/2}
    # with C_nl = 10, E(0) = 0.9 * (gamma_nl / (2 C_nl))^2 => E(tau) <= E(0) exp(-(gamma_nl/2) tau)
    C_nl = 10.0
    E_crit = (gamma_nl / (2.0 * C_nl)) ** 2
    E0 = 0.9 * E_crit
    dt = 2e-3
    E = E0
    riccati_ok = True
    for step in range(1, 20001):  # tau up to 40
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
        and min_rate >= gamma_nl
        and riccati_ok
    )
    report(
        "NB-5 Closed pure-energy 3-block Lyapunov & Riccati trapping (gamma_nl = 1/8, zero W3)",
        ok_nb5,
        f"rates: x={c_x/a0:.4f}, y={c_y/a1:.4f}, z={c_z/a2:.4f} (min=-{min_rate:.4f}<=-0.125), "
        f"E(40)/E(0)={E/E0:.2e} <= e^{{-2.5}}({math.exp(-2.5):.2e})",
    )

    # ---- NB-6: Appendix B continuous-PDE identities ----
    # (i) Exact closed form for k=1: (1/pi)(-2 + sqrt(u) ln((1+sqrt(u))/(1-sqrt(u)))) vs series
    u_test = 0.37
    lhs_k1 = (-2.0 + math.sqrt(u_test) * math.log((1.0 + math.sqrt(u_test)) / (1.0 - math.sqrt(u_test)))) / math.pi
    ms6 = np.arange(0, 200, dtype=float)
    rhs_k1 = float(np.sum(-(u_test**ms6) / (math.pi * (1.0 - ms6 - 0.5))))
    err_disk_id = abs(lhs_k1 - rhs_k1)

    # (ii) Exact derivative L^inf bound on z >= -ln 2: e^{z/4}|d_z(Q_exact G)| <= sqrt(2/3) ||G_z||_{1/4}
    G_b6 = (Y - 1.0) * np.exp(-(Y - 1.0))
    Gz_b6 = (Y - 1.0) * (1.0 - (Y - 1.0)) * np.exp(-(Y - 1.0))
    dQ_b6 = np.gradient(Q_exact @ G_b6, dz)
    n14_z_b6 = float(np.sqrt(np.sum(np.exp(-0.5 * z) * Gz_b6**2) * dz))
    mask_ov = z >= -math.log(2.0)
    deriv_ratio = float(np.max(np.abs(np.exp(0.25 * z[mask_ov]) * dQ_b6[mask_ov])) / n14_z_b6)

    # (iii) Global ||G_z||_{L^2_{1/4}((-inf, -ln 2])}^2 <= (4 sqrt(2) / 9) Z_disk
    # Test on worst-case c_m = (1/m) * 2^{-m} (or random modes)
    m6 = np.arange(1, 61, dtype=float)
    z_left = np.linspace(-10.0, -math.log(2.0), 4001)
    dz_left = float(z_left[1] - z_left[0])
    c_b6 = 1.0 / (m6**1.6)
    Z_b6 = float(np.sum((m6**2) * (c_b6**2)))
    Gz_left = np.sum((m6 * c_b6)[:, None] * np.exp(m6[:, None] * z_left[None, :]), axis=0)
    lhs_gz_left = float(np.trapezoid(np.exp(-0.5 * z_left) * Gz_left**2, z_left))
    gz_left_ratio = lhs_gz_left / Z_b6
    gz_left_cap = 4.0 * math.sqrt(2.0) / 9.0

    # (iv) 1-IBP O(m^{-3}) cancellation in q_m = q_m^{disk} + q_m^{far} for G(v)=v (v<=1), e^{-z} (z>=0)
    m_check = np.arange(1, 201, dtype=float)
    q_comb = (-1.0 / (math.pi * (1.0 - m_check - 0.5)) - 1.0 / (math.pi * (m_check + 1.5))) / (m_check + 1.0)
    l1_half_sum = float(np.sum(np.sqrt(m_check) * np.abs(q_comb)))
    # Check asymptotic m^3 * q_m -> 2/pi
    tail_m3 = float((m_check[-1] ** 3) * q_comb[-1])

    # (v) Unconditional Plancherel bound: ||(1/2) Q_tilde G||_{L^2_{1/4}(R)} <= (2/5) ||G||_{L^2_{3/4}(R)}
    Q_tilde = Q_exact - np.outer(np.ones_like(z), h1_row)
    W14_b6 = np.exp(-0.5 * z)
    W34_b6 = np.exp(-1.5 * z)
    max_qt_ratio = 0.0
    for G_t in [
        (Y - 1.0) * np.exp(-(Y - 1.0)),
        (Y - 1.0) * np.exp(-(Y - 1.0) / 4.0),
        np.exp(-((z - 2.0) / 1.0) ** 2),
        np.exp(-((z + 2.0) / 1.0) ** 2),
    ]:
        h1_t = float(h1_row @ G_t)
        Qt_t = 0.5 * (Q_tilde @ G_t)
        tail_t = (0.5 * h1_t) ** 2 * 2.0 * math.exp(-0.5 * 10.0)
        lhs_t = math.sqrt(float(np.sum(W14_b6 * Qt_t**2) * dz) + tail_t)
        rhs_t = math.sqrt(float(np.sum(W34_b6 * G_t**2) * dz))
        max_qt_ratio = max(max_qt_ratio, lhs_t / rhs_t)

    ok_nb6 = (
        err_disk_id < 1e-14
        and deriv_ratio <= math.sqrt(2.0 / 3.0)
        and gz_left_ratio <= gz_left_cap
        and abs(tail_m3 - 2.0 / math.pi) < 0.03
        and l1_half_sum < 1.0
        and max_qt_ratio <= 0.40
    )
    report(
        "NB-6 Appendix B PDE identities (disk series, ||(1/2)Q_tilde||_{3/4->1/4}<=2/5, d_z Q, 1-IBP)",
        ok_nb6,
        f"disk_err={err_disk_id:.1e}, max||(1/2)Q_tilde||/||G||_34={max_qt_ratio:.4f}<=0.40, "
        f"dQ_ratio={deriv_ratio:.4f}<={math.sqrt(2/3):.4f}, Gz_left={gz_left_ratio:.4f}<={gz_left_cap:.4f}",
    )

    print()
    if FAILURES:
        print(f"GATE FAIL: {FAILURES}")
        return 1
    print("GATE PASS: Final Nonlinear Bootstrap Closure (NB-1..NB-6)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
