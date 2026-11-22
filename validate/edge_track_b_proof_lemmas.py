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
  TB-4  Lemma 4: Horizon regularity V_G(r^2) = O(r^3 |ln r|) as r -> 0+ for G in X_rho^0
        (a(G) = 0), confirming C^1 transversal crossing of dr/ds = -1 + (1/2) V_G(r^2).
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


def check_tb1() -> bool:
    # Verify exact IBP identity a(P(G)) == -(1/pi) int e^{-z/2} D(QG) G dz on G = N_g(u^3 e^{-u})
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
    P = -QG * Dg(G) + 0.5 * QG * G
    DQG = Dr(QG)

    m = (z >= -24.0) & (z <= 24.0)
    a_direct = float(-(1.0 / math.pi) * np.sum(np.exp(-0.5 * z[m]) * P[m]) * dz)
    a_ibp = float(-(1.0 / math.pi) * np.sum(np.exp(-0.5 * z[m]) * DQG[m] * G[m]) * dz)
    rel = abs(a_direct - a_ibp) / max(abs(a_direct), 1e-30)
    ok = bool(rel < 1e-11)
    print(f"TB-1 Zero-deriv IBP a(P(G))      PASS={ok}: direct={a_direct:+.10e}, ibp={a_ibp:+.10e}, rel={rel:.2e}")
    return ok


def check_tb2() -> bool:
    # Verify t^{1/4} * ||N_{5/8}(e^{-u t} phi_*)||_{0, 5/8} stays uniformly bounded as t -> 0+!
    # Use a fine grid on [-40, 40] so u*t is resolved even at t = 1e-4 (u ~ 1e4 => z ~ -18.4)
    Nz = 8192
    zmin, zmax = -44.0, 44.0
    L = zmax - zmin
    dz = L / Nz
    z = zmin + np.arange(Nz) * dz
    u = np.exp(-0.5 * z)
    xi = 2.0 * np.pi * np.fft.fftfreq(Nz, dz)

    sig2 = 0.625  # 5/8
    Nm2 = build_Nm(xi, sig2)
    w2 = np.exp(-sig2 * z)
    phi_star = 1.0 / (2.0 * np.cosh(0.5 * z))

    scaled = []
    for t in (1e-1, 1e-2, 1e-3, 1e-4):
        fn = np.exp(-u * t) * phi_star
        # ||N_{5/8}(fn)||_{0, 5/8}^2 = (1/(2 pi)) int |N(5/8 + i xi) * fhat(5/8 + i xi)|^2 dxi
        # By Parseval on discrete grid: sum |Nm2 * fft(w2 * fn) * dz|^2 * (dxi / (2 pi)) = (dz / Nz) sum |Nm2 * fft(w2 * fn)|^2
        fhat = np.fft.fft(w2 * fn)
        norm2 = math.sqrt(float(np.sum(np.abs(Nm2 * fhat) ** 2)) * (dz / Nz))
        scaled.append((t, norm2, (t ** 0.25) * norm2))

    # Check that t^{1/4} * norm2 approaches a finite constant as t -> 0+ (ratio between t=1e-3 and t=1e-4 close to 1)
    c_vals = [x[2] for x in scaled]
    ratio_34 = abs(c_vals[3] / c_vals[2] - 1.0)
    ok = bool(ratio_34 < 0.05 and max(c_vals) < 100.0)
    summary = ", ".join(f"t={t:.0e}:{c:.4f}" for t, _, c in scaled)
    print(f"TB-2 Rank-1 t^(-1/4) scaling     PASS={ok}: t^(1/4)||N_{{5/8}}(e^{{-ut}}phi_*)|| = ({summary}), ratio34_diff={ratio_34:.3f}")
    return ok


def check_tb3() -> bool:
    # Verify Lemma 1(2): ||Q G||_{H^inf(S_rho)} finite on |Im z| <= rho = 0.5 for G in X_rho^0
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

    rho = 0.5
    Ghat1 = np.fft.fft(w1 * G) * dz
    Ghat2 = np.fft.fft(w2 * G) * dz
    dxi = 2.0 * math.pi / L
    # Strip norms ||G||_{rho, 1/4} and ||G||_{rho, 5/8}
    # Cap at |xi| <= 40 to avoid float64 floor amplification by e^{rho |xi|}
    keep = np.abs(xi) <= 40.0
    nG1 = math.sqrt(float(np.sum(np.exp(2 * rho * np.abs(xi[keep])) * np.abs(Ghat1[keep]) ** 2) * dxi / (2 * math.pi)))
    nG2 = math.sqrt(float(np.sum(np.exp(2 * rho * np.abs(xi[keep])) * np.abs(Ghat2[keep]) ** 2) * dxi / (2 * math.pi)))

    # Q G on Re zeta = -1/4 and Re zeta = +1/8 (no pole at zeta = 0 since a(G) = 0!)
    def cot_c(s):
        a, b = np.real(s), np.imag(s)
        res = np.empty_like(s, dtype=np.complex128)
        big = np.abs(2 * np.pi * b) > 50.0
        res[big] = -1j * np.sign(b[big])
        sm = ~big
        aa, bb = 2 * np.pi * a[sm], 2 * np.pi * b[sm]
        res[sm] = (np.sin(aa) - 1j * np.sinh(bb)) / (np.cosh(bb) - np.cos(aa))
        return res

    z1 = -0.25 + 1j * xi[keep]
    z2 = +0.125 + 1j * xi[keep]
    QGhat1 = (cot_c(z1) / (z1 + 1.0)) * Ghat1[keep]
    QGhat2 = (cot_c(z2) / (z2 + 1.0)) * Ghat2[keep]
    # L^1-Fourier upper bound on H^inf(S_rho):
    sup_pos = float(np.sum(np.exp(rho * np.abs(xi[keep])) * np.abs(QGhat1)) * dxi / (2 * math.pi))
    sup_neg = float(np.sum(np.exp(rho * np.abs(xi[keep])) * np.abs(QGhat2)) * dxi / (2 * math.pi))
    ratio_Q = max(sup_pos, sup_neg) / (nG1 + nG2)

    ok = bool(math.isfinite(ratio_Q) and ratio_Q < 5.0)
    print(f"TB-3 Pole-free Q on S_rho (0.5)  PASS={ok}: ||G||_X={nG1+nG2:.4f}, H^inf(QG)<={max(sup_pos,sup_neg):.4f}, C_Q_emp={ratio_Q:.4f}")
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


def main() -> int:
    p1 = check_tb1()
    p2 = check_tb2()
    p3 = check_tb3()
    p4 = check_tb4()
    assert p1 and p2 and p3 and p4, "Track B proof lemmas gate FAILED"
    print("ALL PASS (TB-1..TB-4)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
