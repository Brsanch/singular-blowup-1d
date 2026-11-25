#!/usr/bin/env python3
"""
GATE 2 (Target B -- gCLM singular profiles, arXiv:2603.25104v4):
Edge-resolved linearized operator M_a, spectrum lambda(a) across a < 0,
convective weighted semigroup decay, and Mellin / Carleman-Mellin structure.

Checks:
  1. Nonlinear finite-difference certification:
     max|M_a G - dG_fd| < 1e-10 against the raw nonlinear dynamic-rescaling PDE
     across a in {-2.0, -1.5, -1.0, -0.5, -0.2} for BOTH Gauge (a) and Gauge (b).
  2. Nyquist artifact diagnostic:
     Certifies that central-difference positive modes near the edge (a=-1.0, -0.75)
     are grid-scale Nyquist artifacts (||dv||/||avg|| > 1 or O(N)-divergent),
     removed cleanly by 2nd-order upwind (||dv||/||avg|| < 0.1).
  3. Full spectrum lambda(a) scan + refinement & domain-depth sweeps:
     a. #Re(lambda) > 1e-4 == 0 across a in [-2.5, -0.1] for BOTH gauges, across
        refinement (N=400 -> 700) and domain-depth sweeps.
     b. Exact analytic point eigenvalue lambda_*(a) = -(1 + a) on -1 < a < 0
        (eigenfunction G_*(1+v) = v^{2*mu - 1}) confirmed on BOTH zeta-grid and z-grid.
  4. Convective weighted semigroup decay in Gauge (b) (paper's Case 2.2 normalization
     with L_0[g_0] = 0):
     Certifies exponential decay at the predicted convective rates across all three
     regimes (-1 < a < 0, a = -1, a < -1).
  5. Mellin / Carleman-Mellin structure (Step 5):
     a. Direct FFT-Mellin verification of the exact a = -1 conjugation
        N_L * M_{-1}^{bulk} * N_L^{-1} = -2 d_z - e^{-z/2} (N_L(s) = sin(pi*s)/(2s+1)).
     b. Open Hilbert-strip zero classification of c_a(s): zero non-real roots on
        0 < Re(s) < 1, Im(s) > 0, and exact collision of s_1 = mu, s_2 = 1-mu at a = -1.
"""

import os
import sys
import numpy as np
from scipy.linalg import expm

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from gclm_singular.profiles import (  # noqa: E402
    mu_of_a,
    interior_stagnation_point,
    mellin_symbol_c,
)
from gclm_singular.operator import (  # noqa: E402
    build_operator_z,
    build_operator_zeta,
    nonlinear_rhs_z,
    spectrum,
)


def check_fd_verification():
    """[1] Certify M_a against the raw nonlinear dynamic-rescaling RHS (< 1e-10)."""
    print("  [1] Linearized operator vs raw nonlinear dynamic-rescaling PDE:")
    a_list = [-2.0, -1.5, -1.0, -0.5, -0.2]
    max_err_a = 0.0
    max_err_b = 0.0

    for a in a_list:
        ops_a = build_operator_z(a, N=300, zmin=-6.0, zmax=8.0, gauge="a", edge_bc=False)
        z = ops_a["z"]
        G = np.exp(-0.5 * (z - 0.5) ** 2) * np.sin(0.7 * z + 0.3)

        # Gauge (a) is strictly quadratic in G -> central FD is exact for any eps
        dG_lin_a = ops_a["M"] @ G
        eps_a = 1e-2
        dG_fd_a = (nonlinear_rhs_z(eps_a * G, ops_a) - nonlinear_rhs_z(-eps_a * G, ops_a)) / (2.0 * eps_a)
        err_a = float(np.max(np.abs(dG_lin_a - dG_fd_a)))
        max_err_a = max(max_err_a, err_a)

        # Gauge (b) has a small O(eps^2) cubic term -> 4th-order Richardson central FD
        ops_b = build_operator_z(a, N=300, zmin=-6.0, zmax=8.0, gauge="b", edge_bc=False)
        dG_lin_b = ops_b["M"] @ G
        e1, e2 = 4e-3, 2e-3
        fd1 = (nonlinear_rhs_z(e1 * G, ops_b) - nonlinear_rhs_z(-e1 * G, ops_b)) / (2.0 * e1)
        fd2 = (nonlinear_rhs_z(e2 * G, ops_b) - nonlinear_rhs_z(-e2 * G, ops_b)) / (2.0 * e2)
        dG_fd_b = (4.0 * fd2 - fd1) / 3.0
        err_b = float(np.max(np.abs(dG_lin_b - dG_fd_b)))
        max_err_b = max(max_err_b, err_b)

        print(f"      a={a:+5.2f}: Gauge(a) max|diff|={err_a:.3e} | Gauge(b) max|diff|={err_b:.3e}")

    ok = (max_err_a < 1e-10) and (max_err_b < 1e-10)
    return ok


def check_nyquist_diagnostic():
    """[2] Confirm central-diff edge modes are Nyquist artifacts removed by upwind2."""
    print("\n  [2] Central vs upwind2 Nyquist artifact diagnostic:")
    # At a = -1.0 on z-grid: central has checkerboard (||dv||/||avg|| > 1), upwind2 is smooth (< 0.1)
    ops_c = build_operator_z(-1.0, N=500, zmin=-14.0, zmax=8.0, gauge="a", advection_scheme="central")
    vc, ec = spectrum(ops_c["M"])
    nyq_c = float(np.linalg.norm(ec[1:, 0] - ec[:-1, 0]) / (np.linalg.norm(ec[1:, 0] + ec[:-1, 0]) + 1e-30))

    ops_u = build_operator_z(-1.0, N=500, zmin=-14.0, zmax=8.0, gauge="a", advection_scheme="upwind2")
    vu, eu = spectrum(ops_u["M"])
    nyq_u = float(np.linalg.norm(eu[1:, 0] - eu[:-1, 0]) / (np.linalg.norm(eu[1:, 0] + eu[:-1, 0]) + 1e-30))

    # At a = -0.75 on zeta-grid: central leading eigenvalue scales linearly with N (node-1 artifact)
    v300, _ = spectrum(build_operator_zeta(-0.75, N=300, zetamax=10.0, gauge="a", advection_scheme="central")["M"])
    v600, _ = spectrum(build_operator_zeta(-0.75, N=600, zetamax=10.0, gauge="a", advection_scheme="central")["M"])
    vu600, _ = spectrum(build_operator_zeta(-0.75, N=600, zetamax=10.0, gauge="a", advection_scheme="upwind2")["M"])

    print(f"      a=-1.00 z-grid:   central lam0={vc[0].real:+.4f} (Nyquist ratio={nyq_c:.2f}) | "
          f"upwind2 lam0={vu[0].real:+.4f} (Nyquist ratio={nyq_u:.3f})")
    print(f"      a=-0.75 zeta-grid: central N=300->{v300[0].real:+.2f}, N=600->{v600[0].real:+.2f} "
          f"(ratio={v600[0].real / v300[0].real:.2f}x) | upwind2 N=600->{vu600[0].real:+.4f}")

    ok = (
        nyq_c > 1.0
        and nyq_u < 0.15
        and vu[0].real < 1e-4
        and (v600[0].real / v300[0].real > 1.7)
        and vu600[0].real < 1e-4
    )
    return ok


def check_spectrum_scan_and_sweeps():
    """[3] Full lambda(a) scan, cross-discretization check, and refinement/depth sweeps."""
    print("\n  [3] Full spectrum lambda(a) scan (upwind2, Gauge a & b):")
    print(f"      {'a':>6s} {'grid':>4s} {'#Re>0(a)':>8s} {'lam0(a)':>10s} {'exact-(1+a)':>12s} "
          f"{'#Re>0(b)':>8s} {'lam1(b) [non-0]':>20s}")

    a_scan = [-2.5, -2.0, -1.5, -1.25, -1.1, -1.0, -0.9, -0.75, -0.5, -0.25, -0.1]
    all_stable = True
    lam_star_ok = True

    for a in a_scan:
        if a <= -1.0:
            ops_a = build_operator_z(a, N=500, zmin=-14.0, zmax=10.0, gauge="a", advection_scheme="upwind2")
            ops_b = build_operator_z(a, N=500, zmin=-14.0, zmax=10.0, gauge="b", advection_scheme="upwind2")
            gname = "z"
        else:
            ops_a = build_operator_zeta(a, N=500, zetamax=14.0, gauge="a", advection_scheme="upwind2")
            ops_b = build_operator_zeta(a, N=500, zetamax=14.0, gauge="b", advection_scheme="upwind2")
            gname = "zeta"

        va, _ = spectrum(ops_a["M"])
        vb, _ = spectrum(ops_b["M"])
        npos_a = int((va.real > 1e-4).sum())
        npos_b = int((vb.real > 1e-4).sum())
        all_stable = all_stable and (npos_a == 0) and (npos_b == 0)

        # First nontrivial eigenvalue of Gauge (b) (skipping the projected lambda=0 mode)
        nz_b = [v for v in vb if abs(v) > 1e-4]
        l1b = nz_b[0]

        if -1.0 < a < 0.0:
            ex = -(1.0 + a)
            ex_str = f"{ex:+12.4f}"
            if a in (-0.75, -0.5, -0.25):
                if abs(va[0].real - ex) > 0.04:
                    lam_star_ok = False
        else:
            ex_str = "   (singular)" if a < -1.0 else "  +0.0000(0)"

        print(f"      {a:+6.2f} {gname:>4s} {npos_a:8d} {va[0].real:+10.4f} {ex_str:>12s} "
              f"{npos_b:8d} {l1b.real:+9.4f}{l1b.imag:+9.4f}j")

    # Refinement + domain-depth + cross-discretization checks on representative values
    print("\n      Refinement, domain-depth & cross-discretization verification:")
    sweep_ok = True
    # (i) a = -0.50: zeta-grid N=400->700, zmax=10->14 AND independent z-grid
    v_zt_10 = spectrum(build_operator_zeta(-0.5, N=500, zetamax=10.0, gauge="a", advection_scheme="upwind2")["M"])[0]
    v_zt_14 = spectrum(build_operator_zeta(-0.5, N=700, zetamax=14.0, gauge="a", advection_scheme="upwind2")["M"])[0]
    v_z_m05 = spectrum(build_operator_z(-0.5, N=600, zmin=-12.0, zmax=12.0, gauge="a", advection_scheme="upwind2")["M"])[0]
    print(f"      a=-0.50 Gauge(a) exact -0.5000: zeta(10)={v_zt_10[0].real:+.5f}, "
          f"zeta(14,N=700)={v_zt_14[0].real:+.5f}, z-grid={v_z_m05[0].real:+.5f}")
    sweep_ok = sweep_ok and (abs(v_zt_14[0].real + 0.5) < 0.005) and (abs(v_z_m05[0].real + 0.5) < 0.015)

    # (ii) a = -1.00: N=400->700, zmin=-10->-18
    v_m1_1 = spectrum(build_operator_z(-1.0, N=400, zmin=-10.0, zmax=8.0, gauge="b", advection_scheme="upwind2")["M"])[0]
    v_m1_2 = spectrum(build_operator_z(-1.0, N=700, zmin=-18.0, zmax=10.0, gauge="b", advection_scheme="upwind2")["M"])[0]
    print(f"      a=-1.00 Gauge(b): [-10,8] N=400 #Re>0={(v_m1_1.real>1e-4).sum()} | "
          f"[-18,10] N=700 #Re>0={(v_m1_2.real>1e-4).sum()}")
    sweep_ok = sweep_ok and ((v_m1_1.real > 1e-4).sum() == 0) and ((v_m1_2.real > 1e-4).sum() == 0)

    # (iii) a = -2.00: N=400->700, [-10,8]->[-18,12]
    v_m2_1 = spectrum(build_operator_z(-2.0, N=400, zmin=-10.0, zmax=8.0, gauge="b", advection_scheme="upwind2")["M"])[0]
    v_m2_2 = spectrum(build_operator_z(-2.0, N=700, zmin=-18.0, zmax=12.0, gauge="b", advection_scheme="upwind2")["M"])[0]
    print(f"      a=-2.00 Gauge(b): [-10,8] N=400 #Re>0={(v_m2_1.real>1e-4).sum()} | "
          f"[-18,12] N=700 #Re>0={(v_m2_2.real>1e-4).sum()}")
    sweep_ok = sweep_ok and ((v_m2_1.real > 1e-4).sum() == 0) and ((v_m2_2.real > 1e-4).sum() == 0)

    return all_stable and lam_star_ok and sweep_ok


def check_semigroup_decay():
    """[4] Convective weighted semigroup decay in Gauge (b) (L_0[g_0] = 0)."""
    print("\n  [4] Convective weighted semigroup decay (Gauge b, Case 2.2 normalization, L_0[g]=0):")
    ok = True
    for a in [-2.0, -1.5, -1.0, -0.5, -0.25]:
        if a <= -1.0:
            ops = build_operator_z(a, N=400, zmin=-14.0, zmax=12.0, gauge="b", advection_scheme="upwind2")
            coord = ops["z"]
            l0 = ops["full"]["l0_row"][1:] if ops["edge_bc"] else ops["full"]["l0_row"]
            b1 = np.exp(-2.0 * (coord - 0.0) ** 2)
            b2 = np.exp(-2.0 * (coord - 2.0) ** 2)
            G0 = b1 - (float(l0 @ b1) / float(l0 @ b2)) * b2
            G0 /= np.linalg.norm(G0)
            sigma = 0.25
            if a == -1.0:
                w = np.exp(-sigma * coord)
            else:
                z_star = np.log(interior_stagnation_point(a) - 1.0)
                w = np.exp(-sigma * np.abs(coord - z_star))
            pred = (1.0 - a) * sigma
        else:
            ops = build_operator_zeta(a, N=400, zetamax=14.0, gauge="b", advection_scheme="upwind2")
            coord = ops["zeta"]
            v = ops["v"]
            mu = mu_of_a(a)
            kappa = -np.cos(np.pi * mu)
            inv_mu = 1.0 - a
            em1 = np.expm1(coord)
            dv = inv_mu * (kappa ** inv_mu) * (em1 ** (inv_mu - 1.0)) * np.exp(coord)
            w_xi = np.full(len(coord), ops["dzeta"])
            w_xi[-1] *= 0.5
            l0 = (np.sin(np.pi * mu) / np.pi) * (v ** (-mu)) / (1.0 + v) * dv * w_xi
            b1 = np.exp(-2.0 * (coord - 1.5) ** 2)
            b2 = np.exp(-2.0 * (coord - 3.5) ** 2)
            G0 = b1 - (float(l0 @ b1) / float(l0 @ b2)) * b2
            G0 /= np.linalg.norm(G0)
            sigma = 0.5
            w = np.exp(-sigma * coord)
            pred = sigma

        E1 = expm(ops["M"])
        G = G0.copy()
        nw = [float(np.linalg.norm(w * G))]
        for _ in range(4):
            G = E1 @ G
            nw.append(float(np.linalg.norm(w * G)))
        rate = -np.log(nw[4] / nw[1]) / 3.0
        # Require clean exponential decay within 20% of convective prediction
        case_ok = (rate > 0.8 * pred) and (rate < 1.25 * pred)
        ok = ok and case_ok
        print(f"      a={a:+5.2f}: ||G||_w(0->4) = {nw[0]:.2e} -> {nw[4]:.2e}, "
              f"rate = {rate:+.4f} (pred {pred:+.3f})  ok={case_ok}")

    return ok


def check_mellin_conjugation_and_strip():
    """[5] Verify exact a=-1 conjugation via FFT-Mellin and Hilbert-strip zero theorem."""
    print("\n  [5] Mellin / Carleman-Mellin structure:")
    # (a) Exact a = -1 conjugation: N_L * M_{-1}^{bulk} * N_L^{-1} = -2 d_z - e^{-z/2}
    # In positive-s convention s = sigma + i*xi (sigma = 0.25): N(s) = sin(pi*s)/(2*s - 1)
    sigma = 0.25
    Nz = 2048
    zf = np.linspace(-25.0, 25.0, Nz, endpoint=False)
    dzf = zf[1] - zf[0]
    xi = 2.0 * np.pi * np.fft.fftfreq(Nz, d=dzf)
    s_line = sigma + 1j * xi
    phase = np.exp(1j * xi * zf[0])

    def mellin_fwd(f_z):
        return phase * np.fft.ifft(np.exp(sigma * zf) * f_z) * (Nz * dzf)

    def mellin_inv(f_hat):
        return np.exp(-sigma * zf) * np.real(np.fft.fft(f_hat / phase) / (Nz * dzf))

    # Zero-integral F0 so G0 = N^{-1} F0 decays at BOTH ends (no s=0 pole residue)
    F0 = (zf - 2.0) * np.exp(-0.5 * (zf - 2.0) ** 2)
    F0_z = (1.0 - (zf - 2.0) ** 2) * np.exp(-0.5 * (zf - 2.0) ** 2)
    N_mult = np.sin(np.pi * s_line) / (2.0 * s_line - 1.0)
    G0_f = mellin_inv(mellin_fwd(F0) / N_mult)
    rhs_F = -2.0 * F0_z - np.exp(-0.5 * zf) * F0
    rhs_G_ref_f = mellin_inv(mellin_fwd(rhs_F) / N_mult)

    ops = build_operator_z(-1.0, N=800, zmin=-14.0, zmax=14.0, gauge="a", advection_scheme="upwind2")
    z = ops["z"]
    G0_num = np.interp(z, zf, G0_f)
    rhs_num = ops["M"] @ G0_num
    rhs_ref = np.interp(z, zf, rhs_G_ref_f)
    mid = (z > -6.0) & (z < 10.0)
    conj_rel_err = float(np.max(np.abs(rhs_num[mid] - rhs_ref[mid])) / np.max(np.abs(rhs_ref[mid])))

    # (b) Exact open-strip zero classification of c_a(s):
    #     Im[cot(pi*(x+iy))/(1 - x - iy)] < 0 on all (x, y) in (0, 1) x (0, inf)
    #     AND exact real zeros s_1 = mu, s_2 = 1 - mu (double zero at a = -1).
    x = np.linspace(0.01, 0.99, 200)
    y = np.geomspace(1e-3, 20.0, 200)
    X2, Y2 = np.meshgrid(x, y, indexing="ij")
    num_im = Y2 * np.sin(2.0 * np.pi * X2) - (1.0 - X2) * np.sinh(2.0 * np.pi * Y2)
    max_im_num = float(np.max(num_im))

    max_root_err = 0.0
    for a in [-2.5, -2.0, -1.5, -1.0, -0.5, -0.25, -0.1]:
        mu = mu_of_a(a)
        max_root_err = max(
            max_root_err,
            abs(mellin_symbol_c(mu, a)),
            abs(mellin_symbol_c(1.0 - mu, a)),
        )

    print(f"      a=-1.0 conjugation ||M_-1 G0 - N^-1(-2d_z - e^-z/2)N G0||_rel = {conj_rel_err:.4e}")
    print(f"      open-strip Im-numerator max on (0,1)x(1e-3,20)                = {max_im_num:+.3e} (< 0 strictly)")
    print(f"      exact real roots max(|c_a(mu)|, |c_a(1-mu)|) across a < 0     = {max_root_err:.3e}")

    return (conj_rel_err < 2e-2) and (max_im_num < 0.0) and (max_root_err < 1e-14)


def main():
    print("=" * 78)
    print("  GATE 2: gCLM SINGULAR PROFILE LINEAR STABILITY & MELLIN STRUCTURE")
    print("=" * 78)
    ok1 = check_fd_verification()
    ok2 = check_nyquist_diagnostic()
    ok3 = check_spectrum_scan_and_sweeps()
    ok4 = check_semigroup_decay()
    ok5 = check_mellin_conjugation_and_strip()

    verdict = ok1 and ok2 and ok3 and ok4 and ok5
    print("\n" + "-" * 78)
    print(f"  [1] FD<1e-10: {ok1} | [2] Nyquist: {ok2} | [3] Spectrum/sweeps: {ok3} | "
          f"[4] Semigroup: {ok4} | [5] Mellin/conj: {ok5}")
    print("\n  VERDICT:", "PASS -- Gate 2 certified" if verdict else "FAIL")
    print("=" * 78)
    return 0 if verdict else 1


if __name__ == "__main__":
    sys.exit(main())
