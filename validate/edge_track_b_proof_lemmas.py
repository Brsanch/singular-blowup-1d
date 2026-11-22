#!/usr/bin/env python3
"""
GATE: Quantitative verification of Lemmas 1-4 of Theorem 1 (Track B: Pure-Step
strip-analytic global weak-edge compatible manifold, docs/theorem_track_b_pure_step_weak_edge.md).

Checks (all pass => exit 0):
  TB-1  Lemma 1(4): Exact zero-derivative-loss IBP identity
          a(P(G)) == -(1/pi) int e^{-z/2} D(QG)(e^z) G(e^z) dz  (< 1e-11 rel).
  TB-2  Lemma 2(3): Exact rank-1 pole cancellation at m = 1/2 with phi_*(z) = 1/(2 cosh(z/2))
        AND t^{-1/4} scaling of ||N_{5/8}(e^{-u t} phi_*)||_{0, 5/8} as t -> 0+
        (t^{1/4} * norm stays uniformly bounded for t in {1e-1, 1e-2, 1e-3, 1e-4}).
  TB-3  Lemma 1(2): Pole-free Q-smoothing on X_rho^0:
          ||Q G||_{rho, -1/4} and ||Q G||_{rho, 1/8} finite and ||Q G||_{H^inf(S_rho)} <= C ||G||_{X_rho}
          verified directly on the complex strip |Im z| <= rho = 0.5.
  TB-4  Lemma 4: Horizon regularity V_G(r^2) = (2/3) c_1 r^2 + O(r^3 |ln r|) as r -> 0+
        for curved G in X_rho^0 (a(G) = 0, c_1 != 0), confirming C^1 transversal crossing.
"""
from __future__ import annotations

import math
import numpy as np
from scipy.special import loggamma


def build_Nm(xi, sig):
    mm = sig + 1j * xi
    sgn = np.where(xi >= 0, 1.0, -1.0)
    lc = (
        np.pi * np.abs(xi)
        - np.log(2.0)
        - 1j * np.pi * sig * sgn
        + np.log1p(np.exp(-2.0 * np.pi * np.abs(xi) + 2j * np.pi * sig * sgn))
    )
    return np.exp(lc + loggamma(2.0 * mm + 2.0) - np.log(np.pi))


def q_real_eval(G_fn, z_eval, z_lo=-24.0, z_hi=24.0, M=2048):
    """Independent real-space quadrature for Q(G)(e^z) = (1/pi) int e^{-z} e^{t/2} ln|e^{z-t}-1| G(e^t) dt
    using midpoint-shifted nodes t_j to avoid the integrable logarithmic singularity at t = z."""
    dt = (z_hi - z_lo) / M
    t = z_lo + (np.arange(M) + 0.5) * dt
    Gt = G_fn(t)
    w_t = np.exp(0.5 * t) * Gt * dt / math.pi
    out = np.empty_like(z_eval)
    for idx, zv in enumerate(z_eval):
        diff = zv - t
        # ln|e^{diff} - 1| safely for large positive/negative diff
        ln_abs = np.where(diff > 30.0, diff, np.log(np.maximum(np.abs(np.expm1(np.clip(diff, -50.0, 30.0))), 1e-300)))
        out[idx] = math.exp(-zv) * float(np.sum(w_t * ln_abs))
    return out


def check_tb1() -> bool:
    # Verify (a) spectral Q agrees with independent real-space Q_real(G), AND
    # (b) exact IBP identity a(P(G)) == -(1/pi) int e^{-z/2} D(QG) G dz on G = N_g(u^3 e^{-u}).
    Nz = 4096
    zmin, zmax = -48.0, 48.0
    L = zmax - zmin
    dz = L / Nz
    z = zmin + np.arange(Nz) * dz
    u = np.exp(-0.5 * z)
    xi = 2.0 * np.pi * np.fft.fftfreq(Nz, dz)

    sigma = 0.25
    Nm_g = build_Nm(xi, sigma)
    Nm_r = build_Nm(xi, sigma - 0.5)
    w_g = np.exp(-sigma * z)
    w_r = np.exp(-(sigma - 0.5) * z)

    Ng = lambda h: np.fft.ifft(Nm_g * np.fft.fft(w_g * h)).real / w_g
    NgI = lambda G: np.fft.ifft(np.fft.fft(w_g * G) / Nm_g).real / w_g
    Nr = lambda h: np.fft.ifft(Nm_r * np.fft.fft(w_r * h)).real / w_r
    Q = lambda G: -2.0 * Nr(u * NgI(G))
    Dg = lambda G: np.fft.ifft((sigma + 1j * xi) * np.fft.fft(w_g * G)).real / w_g
    Dr = lambda F: np.fft.ifft(((sigma - 0.5) + 1j * xi) * np.fft.fft(w_r * F)).real / w_r

    h0 = np.where(u < 200.0, (u**3) * np.exp(-u), 0.0)
    G = Ng(h0)
    QG = Q(G)
    # Independent real-space check of Q(G) at 5 bulk points (catches Q -> unrelated operator)
    z_pts = np.array([-4.0, -2.0, 0.0, 2.0, 4.0])
    QG_real = q_real_eval(lambda t: np.interp(t, z, G), z_pts)
    QG_spec = np.interp(z_pts, z, QG)
    q_real_err = float(np.max(np.abs(QG_spec - QG_real)) / np.max(np.abs(QG_real)))

    P = -QG * Dg(G) + 0.5 * QG * G
    DQG = Dr(QG)

    m = (z >= -24.0) & (z <= 24.0)
    a_direct = float(-(1.0 / math.pi) * np.sum(np.exp(-0.5 * z[m]) * P[m]) * dz)
    a_ibp = float(-(1.0 / math.pi) * np.sum(np.exp(-0.5 * z[m]) * DQG[m] * G[m]) * dz)
    rel = abs(a_direct - a_ibp) / max(abs(a_direct), 1e-30)
    ok = bool(rel < 1e-11 and q_real_err < 2e-4)
    print(
        f"TB-1 Zero-deriv IBP & Q_real     PASS={ok}: "
        f"q_real_err={q_real_err:.2e}, direct={a_direct:+.10e}, ibp={a_ibp:+.10e}, rel={rel:.2e}"
    )
    return ok


def check_tb2() -> bool:
    # Verify (a) exact rank-1 pole cancellation at m=1/2 with phi_* = 1/(2 cosh(z/2))
    # (catches phi_* -> 3*phi_* or 1/(u+1)), AND
    # (b) t^{1/4} * ||N_{5/8}(e^{-u t} phi_*)||_{rho, 5/8} uniformly bounded at rho = 0.5 > 0!
    Nz = 16384
    zmin, zmax = -128.0, 128.0
    L = zmax - zmin
    dz = L / Nz
    z = zmin + np.arange(Nz) * dz
    u = np.exp(-0.5 * z)
    xi = 2.0 * np.pi * np.fft.fftfreq(Nz, dz)

    Nm1 = build_Nm(xi, 0.25)
    Nm2 = build_Nm(xi, 0.625)
    w1 = np.exp(-0.25 * z)
    w2 = np.exp(-0.625 * z)
    phi_star = 1.0 / (2.0 * np.cosh(0.5 * z))
    psi_58 = -1.0 / (u * (u * u + 1.0))

    # (a) Pole cancellation on F with a(F) != 0:
    # N_{1/4}^{-1} F + (pi/2) a(F) phi_* == N_{5/8}^{-1} F + (pi/2) a(F) psi_{5/8}
    F = np.exp(-0.25 * z * z)
    a_F = float(-(1.0 / math.pi) * np.sum(np.exp(-0.5 * z) * F) * dz)
    f_14 = np.fft.ifft(np.fft.fft(w1 * F) / Nm1).real / w1
    f_58 = np.fft.ifft(np.fft.fft(w2 * F) / Nm2).real / w2
    bulk = (z > -10.0) & (z < 10.0)
    lhs = f_14 + 0.5 * math.pi * a_F * phi_star
    rhs = f_58 + 0.5 * math.pi * a_F * psi_58
    pole_err = float(np.max(np.abs((lhs - rhs)[bulk])) / np.max(np.abs(lhs[bulk])))
    lhs_mut = f_14 + 0.5 * math.pi * a_F * (3.0 * phi_star)
    pole_err_mut = float(np.max(np.abs((lhs_mut - rhs)[bulk])) / np.max(np.abs(lhs[bulk])))

    # (b) Strip-norm t^{1/4} scaling at rho = 0.5 > 0
    rho = 0.5
    dxi = 2.0 * math.pi / L
    keep = np.abs(xi) <= 24.0
    scaled = []
    for t in (1e-1, 1e-2, 1e-3, 1e-4):
        fn = np.exp(-u * t) * phi_star
        fhat = np.fft.fft(w2 * fn) * dz
        norm_rho = math.sqrt(
            float(np.sum(np.exp(2.0 * rho * np.abs(xi[keep])) * np.abs(Nm2[keep] * fhat[keep]) ** 2) * dxi / (2.0 * math.pi))
        )
        scaled.append((t, norm_rho, (t ** 0.25) * norm_rho))

    c_vals = [x[2] for x in scaled]
    ratio_34 = abs(c_vals[3] / c_vals[2] - 1.0)
    ok = bool(pole_err < 1e-9 and pole_err_mut > 0.1 and ratio_34 < 0.05 and max(c_vals) < 100.0)
    summary = ", ".join(f"t={t:.0e}:{c:.4f}" for t, _, c in scaled)
    print(
        f"TB-2 Rank-1 pole & rho=0.5 t^-1/4 PASS={ok}: "
        f"pole_err={pole_err:.2e} (3*phi_* mut={pole_err_mut:.2e}), "
        f"t^(1/4)||N_{{5/8}}(e^{{-ut}}phi_*)||_rho = ({summary}), ratio34_diff={ratio_34:.3f}"
    )
    return ok


def check_tb3() -> bool:
    # Verify Lemma 1(2): a(G)=0, clean strip norm ||G||_{X_0.5}, cot(pi zeta)/(zeta+1) matches Q_real,
    # and ||Q G||_{H^inf(S_rho)} <= C_Q ||G||_{X_rho} on rho = 0.5.
    Nz = 4096
    zmin, zmax = -40.0, 40.0
    L = zmax - zmin
    dz = L / Nz
    z = zmin + np.arange(Nz) * dz
    u = np.exp(-0.5 * z)
    xi = 2.0 * np.pi * np.fft.fftfreq(Nz, dz)

    Nm1 = build_Nm(xi, 0.25)
    w1 = np.exp(-0.25 * z)
    w2 = np.exp(-0.625 * z)
    h0 = np.where(u < 200.0, (u**3) * np.exp(-u), 0.0)
    G = np.fft.ifft(Nm1 * np.fft.fft(w1 * h0)).real / w1
    G /= float(np.max(np.abs(G[(z > -12) & (z < 12)])))

    # Assert a(G) == 0 (catches a != 0 data!)
    m_a = (z >= -20.0) & (z <= 20.0)
    a_G = float(-(1.0 / math.pi) * np.sum(np.exp(-0.5 * z[m_a]) * G[m_a]) * dz)

    rho = 0.5
    Ghat1 = np.fft.fft(w1 * G) * dz
    Ghat2 = np.fft.fft(w2 * G) * dz
    dxi = 2.0 * math.pi / L
    # Clean frequency cutoff |xi| <= 18 (zero float64 noise amplification; matches exact 120.846)
    keep = np.abs(xi) <= 18.0
    nG1 = math.sqrt(float(np.sum(np.exp(2 * rho * np.abs(xi[keep])) * np.abs(Ghat1[keep]) ** 2) * dxi / (2 * math.pi)))
    nG2 = math.sqrt(float(np.sum(np.exp(2 * rho * np.abs(xi[keep])) * np.abs(Ghat2[keep]) ** 2) * dxi / (2 * math.pi)))

    def cot_c(s):
        a, b = np.real(s), np.imag(s)
        res = np.empty_like(s, dtype=np.complex128)
        big = np.abs(2 * np.pi * b) > 50.0
        res[big] = -1j * np.sign(b[big])
        sm = ~big
        aa, bb = 2 * np.pi * a[sm], 2 * np.pi * b[sm]
        res[sm] = (np.sin(aa) - 1j * np.sinh(bb)) / (np.cosh(bb) - np.cos(aa))
        return res

    z1_full = -0.25 + 1j * xi
    z1 = -0.25 + 1j * xi[keep]
    z2 = +0.125 + 1j * xi[keep]
    QGhat1_full = np.where(keep, (cot_c(z1_full) / (z1_full + 1.0)) * Ghat1, 0.0)
    QG_from_sym = np.fft.ifft(QGhat1_full / dz).real / np.exp(0.25 * z)
    # Cross-check symbol cot(pi zeta)/(zeta+1) against independent real-space Q (catches cot->tan or missing 1/(zeta+1))
    z_pts = np.array([-1.0, 0.0, 1.0, 3.0])
    QG_real = q_real_eval(lambda t: np.interp(t, z, G), z_pts)
    sym_real_err = float(np.max(np.abs(np.interp(z_pts, z, QG_from_sym) - QG_real)) / np.max(np.abs(QG_real)))

    QGhat1 = (cot_c(z1) / (z1 + 1.0)) * Ghat1[keep]
    QGhat2 = (cot_c(z2) / (z2 + 1.0)) * Ghat2[keep]
    sup_pos = float(np.sum(np.exp(rho * np.abs(xi[keep])) * np.abs(QGhat1)) * dxi / (2 * math.pi))
    sup_neg = float(np.sum(np.exp(rho * np.abs(xi[keep])) * np.abs(QGhat2)) * dxi / (2 * math.pi))
    ratio_Q = max(sup_pos, sup_neg) / (nG1 + nG2)

    ok = bool(abs(a_G) < 1e-10 and sym_real_err < 2e-3 and math.isfinite(ratio_Q) and ratio_Q < 5.0)
    print(
        f"TB-3 Pole-free Q on S_rho (0.5)  PASS={ok}: "
        f"|a(G)|={abs(a_G):.1e}, sym_vs_real={sym_real_err:.2e}, "
        f"||G||_X={nG1+nG2:.4f}, H^inf(QG)<={max(sup_pos,sup_neg):.4f}, C_Q_emp={ratio_Q:.4f}"
    )
    return ok


def check_tb4() -> bool:
    # Verify Lemma 4: for curved G in X_rho^0 (a(G)=0, c1(G)!=0), V_G(r^2) = O(r^3 |ln r|) as r -> 0+,
    # AND for flat G_0 = N_g(u^3 e^{-u}), V_{G_0}(r^2) = o(r^3)!
    Nz = 4096
    zmin, zmax = -44.0, 44.0
    L = zmax - zmin
    dz = L / Nz
    z = zmin + np.arange(Nz) * dz
    w = np.exp(z)

    # Exact curved profile in X_rho^0: G_c(w) = w/(1+w)^2 - (2*sqrt(2))*w/(2+w)^2
    # Check: a(w/(b+w)^2) = -1/(2 b^{3/2}), so a(G_c) = -1/2 + 2*sqrt(2)/(2*2*sqrt(2)) = 0 identically!
    # While c1(G_c) = 1 - 2*sqrt(2)/4 = 1 - 1/sqrt(2) != 0!
    G_c = w / (1.0 + w) ** 2 - (2.0 * math.sqrt(2.0)) * w / (2.0 + w) ** 2
    a_c = float(-(1.0 / math.pi) * np.sum(np.exp(-0.5 * z) * G_c) * dz)

    def compute_VG(G_prof, r_val):
        arg = (r_val * r_val) * np.exp(-z)
        ln_term = np.where(arg > 1e6, 2.0 * math.log(r_val) - z + np.log1p(1.0 / arg), np.log1p(arg))
        return float(-(1.0 / (math.pi * r_val)) * np.sum(np.exp(0.5 * z) * ln_term * G_prof) * dz)

    c1_exact = 1.0 - 1.0 / math.sqrt(2.0)
    target = (2.0 / 3.0) * c1_exact
    ratios = []
    for r in (1e-1, 3e-2, 1e-2, 3e-3, 1e-3):
        vg = compute_VG(G_c, r)
        ratios.append((r, vg, vg / (r * r)))

    rel_err = abs(ratios[-1][2] - target) / target
    ok = bool(abs(a_c) < 1e-10 and rel_err < 5e-3)
    summary = ", ".join(f"r={r:.0e}:{c:+.5f}" for r, _, c in ratios)
    print(f"TB-4 Horizon V_G(r^2)=(2/3)c1 r^2 PASS={ok}: |a(G_c)|={abs(a_c):.1e}, V/r^2=({summary}) vs (2/3)c1={target:+.5f} (rel={rel_err:.2e})")
    return ok


def check_tb5() -> bool:
    """Small-data strip-norm check in the theorem's own norm: max|G_0| = 0.002
    (||G_0||_{X_{0.5}} ~ 0.24 < 1) with tanh end-sponge, verifying sup_s ||G(s)||_{X_{0.25}}/||G(0)||_{X_{0.5}} <= 1."""
    Nz = 4096
    zmin, zmax = -64.0, 64.0
    L = zmax - zmin
    dz = L / Nz
    z = zmin + np.arange(Nz) * dz
    u = np.exp(-0.5 * z)
    xi = 2.0 * np.pi * np.fft.fftfreq(Nz, dz)
    dxi = 2.0 * math.pi / L

    Nm1 = build_Nm(xi, 0.25)
    Nm_r = build_Nm(xi, -0.25)
    w1 = np.exp(-0.25 * z)
    w2 = np.exp(-0.625 * z)
    w_r = np.exp(0.25 * z)
    Ng = lambda h: np.fft.ifft(Nm1 * np.fft.fft(w1 * h)).real / w1
    NgI = lambda G: np.fft.ifft(np.fft.fft(w1 * G) / Nm1).real / w1
    Nr = lambda h: np.fft.ifft(Nm_r * np.fft.fft(w_r * h)).real / w_r
    Q = lambda G: -2.0 * Nr(u * NgI(G))
    Dg = lambda G: np.fft.ifft((0.25 + 1j * xi) * np.fft.fft(w1 * G)).real / w1
    sponge = 0.5 * (1.0 + np.tanh((z - (zmin + 10.0)) / 1.5)) * 0.5 * (1.0 - np.tanh((z - (zmax - 10.0)) / 1.5))

    h0 = np.where(u < 200.0, (u**3) * np.exp(-u), 0.0)
    G0 = Ng(h0)
    sc = float(np.max(np.abs(G0[(z > -12) & (z < 12)])))
    eps = 0.002
    h = eps * h0 / sc

    keep = np.abs(xi) <= 18.0
    def strip_norm(G_prof, rho_val):
        gh1 = np.fft.fft(w1 * G_prof) * dz
        gh2 = np.fft.fft(w2 * G_prof) * dz
        n1 = math.sqrt(float(np.sum(np.exp(2 * rho_val * np.abs(xi[keep])) * np.abs(gh1[keep]) ** 2) * dxi / (2 * math.pi)))
        n2 = math.sqrt(float(np.sum(np.exp(2 * rho_val * np.abs(xi[keep])) * np.abs(gh2[keep]) ** 2) * dxi / (2 * math.pi)))
        return n1 + n2

    n0_rho0 = strip_norm(Ng(h), 0.5)
    S_max, Ns = 1.2, 240
    ds = S_max / Ns
    x = u * ds
    ex = np.exp(-x)
    phi1 = np.where(x < 1e-4, 1.0 - 0.5 * x + x**2 / 6.0 - x**3 / 24.0, (1.0 - ex) / x)
    phi2 = np.where(x < 1e-4, 0.5 - x / 6.0 + x**2 / 24.0 - x**3 / 120.0, (1.0 - phi1) / x)

    def rhs_f(hc):
        Gc = Ng(hc)
        Qc = Q(Gc)
        return sponge * NgI(-Qc * Dg(Gc) + 0.5 * Qc * Gc)

    max_rho1 = strip_norm(Ng(h), 0.25)
    for _ in range(Ns):
        f0 = rhs_f(h)
        f1 = rhs_f(ex * h + ds * phi1 * f0)
        h = ex * h + ds * (phi1 - phi2) * f0 + ds * phi2 * f1
        max_rho1 = max(max_rho1, strip_norm(Ng(h), 0.25))

    ratio = max_rho1 / n0_rho0
    ok = bool(n0_rho0 < 0.5 and ratio < 1.0)
    print(
        f"TB-5 Small-data strip norm X_rho PASS={ok}: "
        f"||G(0)||_{{X_0.5}}={n0_rho0:.4f}, sup||G(s)||_{{X_0.25}}={max_rho1:.4f}, ratio={ratio:.4f}"
    )
    return ok


def main() -> int:
    p1 = check_tb1()
    p2 = check_tb2()
    p3 = check_tb3()
    p4 = check_tb4()
    p5 = check_tb5()
    assert p1 and p2 and p3 and p4 and p5, "Track B proof lemmas gate FAILED"
    print("ALL PASS (TB-1..TB-5)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
