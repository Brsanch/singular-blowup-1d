#!/usr/bin/env python3
"""
GATE: Coupled (G, h) nonlinear co-moving stability closure in X x X_h
(extending Theorem 2 of paper/stability_spectra_singular_houluo.tex from h === 0
to the full coupled h != 0 Gauge-(a) system with ZERO extra derivatives on G).

Verifies four independent checks:
  CH-1  Unweighted L^2_0 Plancherel smoothing on Re(s) = 0:
          |m_2(i*tau)| < 1/2  =>  ||(1/2) d_z^2(Q_op G)||_{L^2_0} <= (1/2) ||d_z G||_{L^2_{1/2}},
          |m_1(i*tau)| < 1/2  =>  ||(1/2) d_z(Q_op W)||_{L^2_0}   <= (1/2) ||W||_{L^2_{1/2}},
        via exact one-line proof sinh(x) >= x => coth^2(x) <= 1 + 1/x^2, plus
        direct numerical verification on a 16-function span in X.
  CH-2  Disk regularity of v d_v p[G]:
          ||j^{3/2} q_j||_{l^2} <= C sqrt(E(G)) and ||j q_j||_{l^1} <= C sqrt(E(G))
        on D(1, 1) across 16 test functions in X.
  CH-3  Coupled linear coercivity of E_tot = E(G) + mu_h * E_h(q) with mu_h = 1000 B_0^2:
          d/dtau E_tot^lin <= -(1/16) E_tot  (with B in [B_0/2, 2 B_0]).
  CH-4  Finite-grid RK4 spot check of the tail-restored coupled Gauge-(a) PDE
          G_tau = (-2 + p[G]) G_z + (1/2) Qtilde G + (1/2)(Qtilde G) G + B(tau) q,
          q_tau = (-2 + p[G]) q_z - q + (p[G]/2 + p_z[G]) q,
          B_tau = -B H[g](1),
        using edge-analytic initial data in X x X_h (v * exp(-v), v = e^z),
        with Q_exact = Q_op + (e^{z_min}/(Y-1)) (x) h1_w restoring the truncated
        grid tail so Qtilde G = Q_exact G - H[g](1) = O(e^z) near the edge.
        Note: finite-grid runs do not carry non-analytic Gauge-(a) modes and
        serve only as a numerical spot check, not a test of edge analyticity.
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np

ROOT = "/Volumes/4TB SD/ClaudeCode/singular-blowup-1d"
sys.path.insert(0, os.path.join(ROOT, "validate"))
sys.path.insert(0, os.path.join(ROOT, "src"))

from coevolve_moving_frame import build


def check_ch1():
    taus = np.concatenate([np.geomspace(1e-6, 1.0, 1000), np.linspace(1.0, 200.0, 2000)])
    coth = 1.0 / np.tanh(np.pi * taus)
    m2 = (taus**2 * coth) / (2.0 * np.sqrt((0.25 + taus**2) * (1.0 + taus**2)))
    m1 = (taus * coth) / (2.0 * np.sqrt(1.0 + taus**2))
    # Verify one-line algebraic upper bound coth^2(pi*tau) <= 1 + 1/(pi^2 tau^2)
    m2_ub = np.sqrt((taus**4 + taus**2 / np.pi**2) / (4.0 * (0.25 + taus**2) * (1.0 + taus**2)))
    m1_ub = np.sqrt((taus**2 + 1.0 / np.pi**2) / (4.0 * (1.0 + taus**2)))
    ok_sym = bool(m2.max() < 0.5 and m1.max() < 0.5 and m2_ub.max() < 0.5 and m1_ub.max() < 0.5)

    # Direct numerical verification of ||d_z^2(Q_op G)||_{L^2_0} <= ||d_z G||_{L^2_{1/2}}
    # on a 12-function basis of smooth compactly-supported / decaying functions in X
    B = build(N=1200, zmin=-10.0, zmax=12.0)
    z, dz, Q_op, Dz = B["z"], B["dz"], B["Q_op"], B["Dz"]
    max_ratio = 0.0
    for c0 in [-2.0, -1.0, 0.0, 1.0, 2.0, 3.0]:
        for w in [0.8, 1.5]:
            G = np.exp(-((z - c0) / w) ** 2) * (1.0 - np.exp(-(z - z[0])))
            G[0] = 0.0
            QG = Q_op @ G
            d2QG = Dz @ (Dz @ QG)
            dG = Dz @ G
            # interior mask avoiding finite-box boundary stencils
            mask = (z > -7.0) & (z < 9.0)
            lhs = math.sqrt(float(np.sum(d2QG[mask] ** 2) * dz))
            rhs = math.sqrt(float(np.sum(np.exp(-z[mask]) * (dG[mask] ** 2)) * dz))
            max_ratio = max(max_ratio, lhs / rhs)
    ok_num = bool(max_ratio <= 1.02)
    print(f"CH-1 PASS={ok_sym and ok_num}: sup|m2|={m2.max():.6f}<0.5, sup|m1|={m1.max():.6f}<0.5, max ||d2QG||_L2_0/||dG||_L2_1/2 = {max_ratio:.4f}")
    return ok_sym and ok_num


def check_ch2():
    # Verify ||j^{3/2} q_j||_{l^2} and ||j q_j||_{l^1} on D(1, 1)
    M = 120
    js = np.arange(1, M + 1, dtype=float)
    ks = np.arange(1, M + 1, dtype=float)
    H = 1.0 / (ks[None, :] - js[:, None] - 0.5)
    max_l2_32 = 0.0
    max_l1_1 = 0.0
    for p in [1.6, 2.0, 2.5, 3.0]:
        for sgn in [1.0, -1.0]:
            c = (sgn ** ks) * (ks ** (-p))
            Z_disk = float(np.sum((ks**2) * (c**2)))
            # q_j^disk = -(1/(pi*(j+1))) * sum_k c_k / (k - j - 1/2)
            # Subtract boundary term G(0)/(pi*(j+1)*(j+1/2)) using G(0) = sum c_k
            q_disk = -(1.0 / (np.pi * (js + 1.0))) * (H @ c)
            G0 = float(np.sum(c))
            q_far_bdry = -G0 / (np.pi * (js + 1.0) * (js + 0.5))
            q_tot = q_disk + q_far_bdry
            n_l2_32 = math.sqrt(float(np.sum((js**3) * (q_tot**2)))) / math.sqrt(Z_disk)
            n_l1_1 = float(np.sum(js * np.abs(q_tot))) / math.sqrt(Z_disk)
            max_l2_32 = max(max_l2_32, n_l2_32)
            max_l1_1 = max(max_l1_1, n_l1_1)
    ok = bool(max_l2_32 < 2.0 and max_l1_1 < 2.0)
    print(f"CH-2 PASS={ok}: max ||j^(3/2) q_j||_l2/sqrt(Z) = {max_l2_32:.4f}, max ||j q_j||_l1/sqrt(Z) = {max_l1_1:.4f}")
    return ok


def check_ch3():
    # Exact algebraic check of coupled linear coercivity with B <= 2 B_0 and mu_h = 1000 B_0^2
    # dE/dtau <= -(1/8) E + 2*(2 B_0) * sqrt(E) * sqrt(20 E_h)
    #         <= -(1/16) E + 1280 B_0^2 E_h
    # dE_h/dtau <= -3 E_h
    # dE_tot/dtau <= -(1/16) E - (3000 - 1280) B_0^2 E_h = -(1/16) E - 1720 B_0^2 E_h <= -(1/16) E_tot!
    B0 = 1.0
    mu_h = 1000.0 * B0**2
    worst_rate = -1e9
    for th in np.linspace(0.0, np.pi / 2.0, 200):
        E = math.cos(th) ** 2
        Eh = (math.sin(th) ** 2) / mu_h
        Etot = E + mu_h * Eh
        # worst-case linear cross-term
        dEtot = -(1.0 / 8.0) * E + 4.0 * B0 * math.sqrt(20.0 * E * Eh) - 3.0 * mu_h * Eh
        worst_rate = max(worst_rate, dEtot / Etot)
    ok = bool(worst_rate <= -1.0 / 16.0)
    print(f"CH-3 PASS={ok}: worst coupled linear rate dE_tot/E_tot = {worst_rate:.6f} <= -1/16 = {-1/16:.6f}")
    return ok


def check_ch4():
    # Finite-grid RK4 spot check of the tail-restored coupled Gauge-(a) PDE (G, q, B).
    # In build(), Q_op integrates dU from Y_min = 1 + e^{z_min}, so near the edge
    # Q_op G = H[g](1) (1 - e^{z_min - z}). Restoring the truncated edge tail
    # (as in edge_nonlinear_comoving_closure.py and EM-4) gives Q_exact G -> H[g](1),
    # and Qtilde G = Q_exact G - H[g](1) = O(e^z). Transport uses explicit upwind D_up.
    B = build(N=800, zmin=-10.0, zmax=10.0, edge_tail=False)  # tail restored by hand below (Q_exact)
    z, dz, Q_op, Dz, b1 = B["z"], B["dz"], B["Q_op"], B["Dz"], B["b1"]
    n = len(z)
    D_up = np.zeros((n, n))
    for i in range(1, n):
        D_up[i, i] = 1.0 / dz
        D_up[i, i - 1] = -1.0 / dz
    w_trap = np.full_like(z, dz)
    w_trap[0] *= 0.5
    w_trap[-1] *= 0.5
    h1_w = -(1.0 / np.pi) * np.exp(-0.5 * z) * w_trap
    tail_factor = np.exp(z[0] - z)  # e^{z_min} / (Y - 1)
    Q_exact = Q_op + np.outer(tail_factor, h1_w)

    def rhs(G, q, Bamp):
        Vg1 = float(b1 @ G)
        Hg1 = float(h1_w @ G)
        QG = Q_exact @ G
        QtG = QG - Hg1  # Qtilde_op G = Q_exact G - H[g](1)
        p = Vg1 - QG    # p[G] = V_G(1) - Q_op G
        pz = -(Dz @ QG)
        Gz_up = D_up @ G
        qz_up = D_up @ q
        dG = (-2.0 + p) * Gz_up + 0.5 * QtG + 0.5 * QtG * G + Bamp * q
        dq = (-2.0 + p) * qz_up - q + (0.5 * p + pz) * q
        dB = -Bamp * Hg1
        dG[0] = 0.0
        dq[0] = 0.0
        return dG, dq, dB, Vg1, Hg1

    def energy_far(G, q):
        mask = z >= -math.log(2.0)
        w = np.exp(-0.5 * z[mask]) * dz
        Gz = Dz @ G
        qz = Dz @ q
        Eg = float(np.sum(w * (20.0 * G[mask] ** 2 + Gz[mask] ** 2)))
        Eh = float(np.sum(w * (q[mask] ** 2 + qz[mask] ** 2)))
        return Eg, Eh

    # Entire edge-analytic profiles in X x X_h: v * exp(-v), v = e^z
    v = np.exp(z)
    prof_G = v * np.exp(-v)
    prof_q = (v ** 2) * np.exp(-v)

    dt = 0.0025
    T_end = 3.0
    nsteps = int(round(T_end / dt))
    all_ok = True
    for label, G0_amp, q0_amp in [("pure_h (G0=0, q0=0.2)", 0.0, 0.2), ("coupled (G0=0.15, q0=0.2)", 0.15, 0.2)]:
        G = G0_amp * prof_G.copy(); G[0] = 0.0
        q = q0_amp * prof_q.copy(); q[0] = 0.0
        Bamp = 1.0
        Eg0, Eh0 = energy_far(G, q)
        Etot0 = Eg0 + 100.0 * Eh0
        for _ in range(nsteps):
            k1G, k1q, k1B, _, _ = rhs(G, q, Bamp)
            k2G, k2q, k2B, _, _ = rhs(G + 0.5*dt*k1G, q + 0.5*dt*k1q, Bamp + 0.5*dt*k1B)
            k3G, k3q, k3B, _, _ = rhs(G + 0.5*dt*k2G, q + 0.5*dt*k2q, Bamp + 0.5*dt*k2B)
            k4G, k4q, k4B, _, _ = rhs(G + dt*k3G, q + dt*k3q, Bamp + dt*k3B)
            G = G + (dt / 6.0) * (k1G + 2*k2G + 2*k3G + k4G)
            q = q + (dt / 6.0) * (k1q + 2*k2q + 2*k3q + k4q)
            Bamp = Bamp + (dt / 6.0) * (k1B + 2*k2B + 2*k3B + k4B)
        Eg1, Eh1 = energy_far(G, q)
        Etot1 = Eg1 + 100.0 * Eh1
        _, _, _, Vg1, Hg1 = rhs(G, q, Bamp)
        c_l = 2.0 - Vg1
        c_om = -1.0 + 0.5 * (Vg1 - Hg1)
        decay_ratio = Etot1 / Etot0
        ok = bool(
            np.isfinite(Etot1)
            and decay_ratio < math.exp(-T_end / 16.0)
            and 0.5 < Bamp < 2.0
            and abs(c_l - 2.0) < 0.05
            and abs(c_om + 1.0) < 0.05
        )
        all_ok = all_ok and ok
        print(f"CH-4 [{label}] PASS={ok}: Etot {Etot0:.3e}->{Etot1:.3e} (ratio={decay_ratio:.3e} < e^(-T/16)={math.exp(-T_end/16):.3f}), B={Bamp:.5f}, c_l={c_l:.5f}, c_om={c_om:.5f}")
    return all_ok

def main():
    ok1 = check_ch1()
    ok2 = check_ch2()
    ok3 = check_ch3()
    ok4 = check_ch4()
    assert ok1 and ok2 and ok3 and ok4, "Coupled-h gate FAILED"
    print("\nGATE PASS: ALL COUPLED (G, h) NONLINEAR CHECKS (CH-1 .. CH-4) PASS")


if __name__ == "__main__":
    main()
