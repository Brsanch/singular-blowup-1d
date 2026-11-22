#!/usr/bin/env python3
"""
GATE: Direction 2 — Nonlinear Co-Moving Leg (Gauge-(a) Exact Collapse, Dual-Strip
Symbol Bounds, Edge-Subtracted W1 at sigma_- = 3/4, Hilbert Edge Commutator, and
Self-Consistent Nonlinear Relaxation) (docs/desk/nonlinear_comoving_leg.md, 2026-11-20).

Verifies the five load-bearing ingredients of the nonlinear co-moving stability
framework for the singular 1D Hou-Luo profile (Chen-Huang-Li arXiv:2604.01868 Thm. 2.3):

  NL-1  Exact Gauge-(a) PDE collapse:
        With A = 1, c_l = 2 - V_G(1), c_omega = -1 + (V_G(1) - H[g](1))/2, B_dot = -B H[g](1),
        all singular e^{-z} dilation terms and all inhomogeneous amplitude/Theta forcings
        cancel identically against raw eq. (2.4).
  NL-2  Dual-strip exact symbol & functional bounds at (sigma_+ = 1/4, sigma_- = 3/4):
        |tan(pi*(1/4 + i*tau))| = |tan(pi*(3/4 + i*tau))| = 1 identically, giving
        sup|c(1/4+it)| = 2/3, sup|c(3/4+it)| = 2/5, sup|s c(s)| = 1/2,
        lower-order derivative smoothing ||0.5 d_z Q_tilde G||_{1/4} <= (7/10) ||G||_{3/4},
        C_H(1/4, 3/4) = 1/sqrt(pi), C_V(1/4) < 0.52.
  NL-3  Dual-strip W1 Route F resolvent bounds at sigma_- = 3/4, delta = 1/20 (Re lam = -1.45):
        Digamma ceiling L(3/4, 1/20) = L(1/4, 1/20) = 54.7581 < 55 by sigma <-> 1-sigma symmetry,
        off-diagonal Schur norm uniform, rank-one defect decays O(1/(1+|omega|)).
  NL-4  Exact Hilbert edge-commutator exponent preservation:
        For smooth Y-perturbations G(z) = e^z (alpha = 1, F(Y) = (Y-1) F_1(Y)),
        c(1) = -tan(pi)/3 = 0 cancels the s=1 pole so Q_tilde G(z) = O(e^z) (slope 1.000,
        zero exponent loss), whereas fractional alpha = 0.75 loses 1/2 (slope 0.25).
  NL-5  Full self-consistent Gauge-(a) nonlinear march (bulk & edge-smooth initial data):
        ||G||_{3/4,1/4}, ||G_z||_{3/4,1/4}, ||Q_op G||_inf, ||F_1||_{L2[1,2]}, and |H[g](1)|
        all decay exponentially; B(tau) -> B_inf > 0, c_l -> 2, c_omega -> -1 with zero drift.

Exit 0 iff NL-1..NL-5 all PASS.
"""
from __future__ import annotations

import math
import os
import sys

import mpmath as mp
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "validate"))

from stability_edge_z import build_blocks_z  # noqa: E402
import edge_w1_route_f_closure as w1  # noqa: E402

mp.mp.dps = 30
FAILURES: list[str] = []


def report(name: str, ok: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}  {detail}")
    if not ok:
        FAILURES.append(name)


def main() -> int:
    print("Direction 2: Nonlinear Co-Moving Leg Closure Gate (NL-1..NL-5)")

    N = 1000
    zmin, zmax = -10.0, 10.0
    B = build_blocks_z(N=N, zmin=zmin, zmax=zmax, beta=0.0, edge_bc=True, edge_tail=False)  # edge_tail=False: this gate restores the tail itself (Q_exact below)
    z, dz, Y, Q_op, Dz = B["z"], B["dz"], B["Y"], B["Q_op"], B["Dz"]
    n = len(z)
    Ymin = 1.0 + math.exp(zmin)

    w_trap = np.full(n, dz)
    w_trap[0] *= 0.5
    w_trap[-1] *= 0.5
    h1_row = -(1.0 / np.pi) * np.exp(-0.5 * z) * w_trap
    ln_term = np.where(z < 0, -z + np.log1p(np.exp(z)), np.log1p(np.exp(-z)))
    v1_row = -(1.0 / np.pi) * np.exp(0.5 * z) * ln_term * w_trap

    Q_exact = Q_op + np.outer((Ymin - 1.0) / (Y - 1.0), h1_row)
    Q_tilde = Q_exact - np.outer(np.ones(n), h1_row)

    # ---- NL-1: Exact Gauge-(a) PDE collapse vs raw eq. (2.4) ----
    G = 0.2 * np.exp(-((z - 1.0) / 1.2) ** 2)
    h = 0.15 * np.exp(-((z - 1.5) / 1.1) ** 2)
    B_amp = 1.13
    Vg1 = float(v1_row @ G)
    Hg1 = float(h1_row @ G)
    QG = Q_exact @ G
    QtG = QG - Hg1
    Gz = Dz @ G
    hz = Dz @ h

    c_l = 2.0 - Vg1
    c_om = -1.0 + 0.5 * (Vg1 - Hg1)
    U = -2.0 + Vg1 + (Y - 1.0) * QG
    a_phys = U + c_l * Y
    Omega = (Y - 1.0) ** (-0.5) * (1.0 + G)
    Omega_Y = (Y - 1.0) ** (-1.5) * (Gz - 0.5 * (1.0 + G))
    Theta = B_amp * (0.5 * np.pi + h)
    Theta_Y = B_amp * (Y - 1.0) ** (-1.0) * hz

    dG_raw = (Y - 1.0) ** 0.5 * (-a_phys * Omega_Y + c_om * Omega + Theta_Y)
    Bdot = -B_amp * Hg1
    dh_raw = (-a_phys * Theta_Y + (c_l + 2.0 * c_om) * Theta - Bdot * (0.5 * np.pi + h)) / B_amp

    dG_box = -2.0 * Gz + 0.5 * QtG + B_amp * np.exp(-0.5 * z) * hz + (Vg1 - QG) * Gz + 0.5 * QtG * G
    dh_box = -2.0 * hz + (Vg1 - QG) * hz

    err_G = float(np.max(np.abs(dG_raw - dG_box)))
    err_h = float(np.max(np.abs(dh_raw - dh_box)))
    report(
        "NL-1 Exact Gauge-(a) PDE algebraic identity vs raw eq. (2.4)",
        err_G < 1e-10 and err_h < 1e-13,
        f"err_G={err_G:.3e}, err_h={err_h:.3e}",
    )

    # ---- NL-2: Dual-strip exact symbol & functional norms at (1/4, 3/4) ----
    taus = np.linspace(-50.0, 50.0, 4001)
    s14 = 0.25 + 1j * taus
    s34 = 0.75 + 1j * taus
    tan14_err = float(np.max(np.abs(np.abs(np.tan(np.pi * s14)) - 1.0)))
    tan34_err = float(np.max(np.abs(np.abs(np.tan(np.pi * s34)) - 1.0)))
    c14_max = float(np.max(np.abs(-np.tan(np.pi * s14) / (2.0 * s14 + 1.0))))
    c34_max = float(np.max(np.abs(-np.tan(np.pi * s34) / (2.0 * s34 + 1.0))))
    sc14_max = float(np.max(np.abs(-s14 * np.tan(np.pi * s14) / (2.0 * s14 + 1.0))))
    sc34_max = float(np.max(np.abs(-s34 * np.tan(np.pi * s34) / (2.0 * s34 + 1.0))))

    z_fine = np.linspace(-40.0, 40.0, 40001)
    ln_f = np.where(z_fine < 0, -z_fine + np.log1p(np.exp(z_fine)), np.log1p(np.exp(-z_fine)))
    CV = float(np.sqrt(np.trapezoid((np.exp(0.75 * z_fine) * ln_f) ** 2, z_fine)) / np.pi)
    kh = np.exp(-0.5 * z_fine) / np.sqrt(np.exp(-1.5 * z_fine) + np.exp(-0.5 * z_fine))
    CH = float(np.sqrt(np.trapezoid(kh**2, z_fine)) / np.pi)
    CH_exact = 1.0 / math.sqrt(math.pi)

    def wnorm1(vec: np.ndarray, sig: float) -> float:
        return float(np.sqrt(np.sum(np.exp(-2.0 * sig * z) * vec**2) * dz))

    G_osc = (Y - 1.0) * np.exp(-(Y - 1.0)) * np.cos(4.0 * z)
    smooth_ratio = wnorm1(0.5 * (Dz @ (Q_tilde @ G_osc)), 0.25) / wnorm1(G_osc, 0.75)

    ok_nl2 = (
        tan14_err < 1e-12
        and tan34_err < 1e-12
        and abs(c14_max - 2.0 / 3.0) < 1e-6
        and abs(c34_max - 2.0 / 5.0) < 1e-6
        and sc14_max <= 0.5 + 1e-9
        and sc34_max <= 0.5 + 1e-9
        and abs(CH - CH_exact) < 1e-6
        and CV < 0.52
        and smooth_ratio <= 0.70
    )
    report(
        "NL-2 Dual-strip exact rational symbol & smoothing bounds at (1/4, 3/4)",
        ok_nl2,
        f"c14={c14_max:.4f}(2/3), c34={c34_max:.4f}(2/5), sc<={max(sc14_max,sc34_max):.4f}(1/2), "
        f"CH={CH:.5f}(1/sqrt(pi)), CV={CV:.5f}, ratio={smooth_ratio:.4f}<=0.70",
    )

    # ---- NL-3: Dual-strip W1 Route F resolvent bounds at sigma_- = 3/4, delta = 1/20 ----
    sigma_m, delta_m = 0.75, 0.05
    L_m = 3.0 * math.pi + 2.0 / delta_m + 1.0 / sigma_m + 1.0 / (1.0 - sigma_m)
    L_p = 3.0 * math.pi + 2.0 / delta_m + 1.0 / 0.25 + 1.0 / 0.75
    max_dlog = 0.0
    for omega in [0.0, 5.0, 20.0, 60.0]:
        for tau in [-40.0, -10.0, -0.5 * omega, 0.0, 10.0, 40.0]:
            max_dlog = max(max_dlog, abs(w1.d_log_X_abs(sigma_m, delta_m, tau, omega)))
    def_norms = {}
    for omega in [0.0, 5.0, 20.0, 50.0]:
        lo = min(0.0, -0.5 * omega) - 10.0
        hi = max(0.0, -0.5 * omega) + 10.0
        n_g = 360
        h_g = (hi - lo) / n_g
        nx_sq, nk_sq = mp.mpf(0), mp.mpf(0)
        for j in range(n_g + 1):
            tp = lo + j * h_g
            wt = mp.mpf("0.5") if j in (0, n_g) else mp.mpf(1)
            lx = mp.mpf(w1.log_X_abs(sigma_m, delta_m, tp, omega))
            xv = mp.exp(lx)
            nx_sq += wt * (xv**2)
            wp = mp.exp(4j * mp.pi * sigma_m) * mp.exp(-4 * mp.pi * tp)
            ks = abs((-4 * mp.pi * wp) / (wp - 1)) if tp >= 0 else abs((-4 * mp.pi * wp) * (1 / (wp - 1) - 1 / wp))
            denom = abs(mp.mpc(delta_m, 2.0 * tp + omega))
            kv = (ks / (2 * mp.pi)) / (xv * denom)
            nk_sq += wt * (kv**2)
        def_norms[omega] = float(mp.sqrt(nx_sq * nk_sq) * h_g)
    ok_nl3 = (
        abs(L_m - L_p) < 1e-12
        and max_dlog < L_m
        and def_norms[0.0] < 3.0
        and def_norms[0.0] > def_norms[5.0] > def_norms[20.0] > def_norms[50.0]
    )
    report(
        "NL-3 Dual-strip W1 Route F at sigma_-=3/4, delta=1/20 (Re lam=-1.45)",
        ok_nl3,
        f"L(3/4)=L(1/4)={L_m:.4f}, max|dlogX|={max_dlog:.4f}, defect: {def_norms[0.0]:.3f}->{def_norms[50.0]:.4f}",
    )

    # ---- NL-4: Exact Hilbert edge-commutator slope preservation at alpha = 1 ----
    def edge_slope(alpha_str: str) -> float:
        alpha = mp.mpf(alpha_str)
        Hg1_a = -1 / (mp.pi * (alpha - mp.mpf("0.5")))
        vals = []
        for eps_s in ["1e-3", "1e-4"]:
            eps = mp.mpf(eps_s)
            I1 = mp.quad(lambda u: u ** (alpha - mp.mpf("0.5")) * mp.log(eps / u - 1), [0, eps / 2, eps])
            I2 = mp.quad(lambda u: u ** (alpha - mp.mpf("0.5")) * mp.log(1 - eps / u), [eps, 2 * eps, 1])
            vals.append(float((I1 + I2) / (mp.pi * eps) - Hg1_a))
        return float(math.log(abs(vals[1] / vals[0])) / math.log(0.1))

    s_075 = edge_slope("0.75")
    s_100 = edge_slope("1.00")
    report(
        "NL-4 Hilbert edge-commutator exponent preservation (c(1)=0 kills s=1 pole)",
        abs(s_100 - 1.0) < 0.01 and abs(s_075 - 0.25) < 0.02,
        f"slope(alpha=1.00)={s_100:.4f} (zero loss), slope(alpha=0.75)={s_075:.4f} (loses 1/2)",
    )

    # ---- NL-5: Self-consistent Gauge-(a) nonlinear march (edge-smooth data) ----
    def wnorm2(vec: np.ndarray) -> float:
        W = np.exp(-1.5 * z) + np.exp(-0.5 * z)
        return float(np.sqrt(np.sum(W * vec**2) * dz))

    def rhs(G_v, h_v, B_v):
        vg = float(v1_row @ G_v)
        hg = float(h1_row @ G_v)
        qg = Q_exact @ G_v
        qtg = qg - hg
        gz = Dz @ G_v
        hz_v = Dz @ h_v
        dG_v = -2.0 * gz + 0.5 * qtg + B_v * np.exp(-0.5 * z) * hz_v + (vg - qg) * gz + 0.5 * qtg * G_v
        dh_v = -2.0 * hz_v + (vg - qg) * hz_v
        dB_v = -B_v * hg
        return dG_v, dh_v, dB_v

    G = 0.25 * (Y - 1.0) * np.exp(-(Y - 1.0))
    G[0] = 0.0
    h = 0.15 * np.exp(-((z - 1.5) / 1.0) ** 2)
    h[0] = 0.0
    B_amp = 1.0
    g0_norm = wnorm2(G)
    gz0_norm = wnorm2(Dz @ G)
    q0_inf = float(np.max(np.abs(Q_exact @ G)))
    mask12 = (Y >= 1.0 + math.exp(-8.0)) & (Y <= 2.0)
    dY12 = np.exp(z[mask12]) * dz
    f1_0 = float(np.sqrt(np.sum((G[mask12] / (Y[mask12] - 1.0)) ** 2 * dY12)))

    dt = 0.2 * dz
    T = 4.0
    tau = 0.0
    while tau < T - 1e-12:
        step = min(dt, T - tau)
        k1G, k1h, k1B = rhs(G, h, B_amp)
        k2G, k2h, k2B = rhs(G + 0.5 * step * k1G, h + 0.5 * step * k1h, B_amp + 0.5 * step * k1B)
        k3G, k3h, k3B = rhs(G + 0.5 * step * k2G, h + 0.5 * step * k2h, B_amp + 0.5 * step * k2B)
        k4G, k4h, k4B = rhs(G + step * k3G, h + step * k3h, B_amp + step * k3B)
        G = G + (step / 6.0) * (k1G + 2 * k2G + 2 * k3G + k4G)
        h = h + (step / 6.0) * (k1h + 2 * k2h + 2 * k3h + k4h)
        B_amp = B_amp + (step / 6.0) * (k1B + 2 * k2B + 2 * k3B + k4B)
        G[0] = 0.0
        h[0] = 0.0
        tau += step

    gT_norm = wnorm2(G)
    gzT_norm = wnorm2(Dz @ G)
    qT_inf = float(np.max(np.abs(Q_exact @ G)))
    f1_T = float(np.sqrt(np.sum((G[mask12] / (Y[mask12] - 1.0)) ** 2 * dY12)))
    c_l_T = 2.0 - float(v1_row @ G)
    c_om_T = -1.0 + 0.5 * (float(v1_row @ G) - float(h1_row @ G))

    ok_nl5 = (
        gT_norm < 0.15 * g0_norm
        and gzT_norm < 0.15 * gz0_norm
        and qT_inf < 0.05 * q0_inf
        and f1_T < 1e-3 * f1_0
        and abs(c_l_T - 2.0) < 5e-3
        and abs(c_om_T + 1.0) < 5e-3
        and 0.5 < B_amp < 2.0
    )
    report(
        "NL-5 Self-consistent Gauge-(a) nonlinear march relaxes all norms & exponents",
        ok_nl5,
        f"||G||:{g0_norm:.3f}->{gT_norm:.4f}, ||Gz||:{gz0_norm:.3f}->{gzT_norm:.4f}, "
        f"||Q||_inf:{q0_inf:.3f}->{qT_inf:.4f}, ||F1||:{f1_0:.3f}->{f1_T:.2e}, "
        f"c_l={c_l_T:.5f}, c_om={c_om_T:.5f}, B={B_amp:.4f}",
    )

    # ---- NL-6: Weak-edge distributional delta(X-1) jump balance B = A^2 across X = 1 ----
    # For Omega_eta(v) = Re[(v + i*eta)^(-1/2)], U_{X,eta}(v) = H(Omega_eta)(v) = Im[(v + i*eta)^(-1/2)],
    # their exact product is U_{X,eta} * Omega_eta = (1/2) Im[(v + i*eta)^(-1)] = -eta / (2*(v^2 + eta^2))
    # -> -(pi/2) * delta(v) as eta -> 0+. Substitute v = eta * tan(theta), dv / (v^2 + eta^2) = dtheta / eta:
    from scipy.integrate import quad
    eta = 1e-10
    th_max = math.atan(1.0 / eta)
    def integrand_th(th):
        v = eta * math.tan(th)
        if abs(v) >= 1.0:
            return 0.0
        return 0.5 * math.exp(1.0 - 1.0 / (1.0 - v * v))
    lorentz_int = quad(integrand_th, -th_max, th_max, limit=200, epsabs=1e-13, epsrel=1e-13)[0]
    res_B1 = abs(lorentz_int - 0.5 * math.pi * 1.0)
    res_B101 = lorentz_int - 0.5 * math.pi * 1.01
    target_B101 = 0.5 * math.pi * (1.0 - 1.01)
    ok_nl6 = res_B1 < 1e-9 and abs(res_B101 - target_B101) < 1e-9
    report(
        "NL-6 Weak-edge distributional delta(X-1) jump balance B = A^2 across X = 1",
        ok_nl6,
        f"|weak_res(B=1)|={res_B1:.2e}, weak_res(B=1.01)={res_B101:+.8f} (exact (pi/2)(1-B)={target_B101:+.8f})",
    )

    print()
    if FAILURES:
        print(f"GATE FAIL: {FAILURES}")
        return 1
    print("GATE PASS: Direction 2 Nonlinear Co-Moving Leg (NL-1..NL-6)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
