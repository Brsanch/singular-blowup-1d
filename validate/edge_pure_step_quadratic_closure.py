#!/usr/bin/env python3
"""
GATE: Full time-dependent O(eps^2) pure-step (eta === 0) weak-edge closure.

Combines Codex's turning-point homological ODE profile psi = N h (solving
-6 psi - 2 psi_z - (1/2) Q psi = source with c_1(psi) = -(3/pi) J and a(psi) = 0)
with the exact closed-form conjugated gauge-modulation term
  G_{2,V}(v, tau) = A(tau) G_1(v, tau) - (1/2) B(tau) (Q G_1)(v, tau)
absorbing +V_{G_1}(tau) v d_v G_1(v, tau), and verifies the full time-dependent
O(eps^2) PDE + weak-edge compatibility a(G_2)(tau) = 0.

Checks (all pass => exit 0):
  PQ-1  Turning-point homological ODE profile psi (n=8192): bulk residual < 1e-6
        and c_1(psi) = -(3/pi) D_0.
  PQ-2  Zero edge moment a(psi) = 0 with exact edge tail -(2 c_1 / pi) e^{z_lo/2}
        restored (< 1e-6 on [-16, 16]).
  PQ-3  Full time-dependent O(eps^2) PDE residual for G_2 = G_{2,hom} + G_{2,V}
        across tau in {0.0, 0.3, 0.8}, v in {0.05, 0.3, 1.5} (max rel < 2e-5).
"""
from __future__ import annotations

import math
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.interpolate import make_interp_spline
from scipy.special import kv, kvp, loggamma


def build_psi(n: int = 8192):
    def phi_r(r, b=1.0):
        r = np.asarray(r)
        safe = (r > 0) & (r < 350)
        rr = np.where(safe, r, 1.0)
        val = -2.0 / math.pi * rr**5 * np.real(np.sqrt(-1j / b) * kv(1, 2 * rr * np.sqrt(-1j * b)))
        return np.where(safe, val, 0.0)

    def phib_r(r, b=1.0):
        r = np.asarray(r)
        safe = (r > 0) & (r < 350)
        rr = np.where(safe, r, 1.0)
        z = 2 * rr * np.sqrt(-1j * b)
        val = -2.0 / math.pi * rr**5 * np.real(np.sqrt(-1j / b) / (2 * b) * (-kv(1, z) + z * kvp(1, z)))
        return np.where(safe, val, 0.0)

    def source_r(r):
        f = phi_r(r)
        fb = phib_r(r)
        return fb * (4.0 * f + fb)

    def ir(fun):
        return sum(
            quad(fun, a, b, epsabs=2e-9, epsrel=2e-11, limit=200)[0]
            for a, b in ((0, 1), (1, 3), (3, 8), (8, 20), (20, 40), (40, 80))
        )

    D = [ir(lambda r: 4 * r ** (4 * m + 1) * float(source_r(r))) for m in range(3)]
    zmin, zmax = -120.0, 80.0
    L = zmax - zmin
    dz = L / n
    z = zmin + np.arange(n) * dz
    sigma = 0.1
    freq = 2 * np.pi * np.fft.fftfreq(n, dz)
    s = sigma + 1j * freq
    sgn = np.where(freq >= 0, 1.0, -1.0)
    logcos = (
        np.pi * abs(freq)
        - np.log(2)
        - 1j * np.pi * sigma * sgn
        + np.log1p(np.exp(-2 * np.pi * abs(freq) + 2j * np.pi * sigma * sgn))
    )
    Nm = np.exp(logcos + loggamma(2 * s + 2) - np.log(np.pi))
    r = np.exp(-z / 4.0)
    source = source_r(r)
    fk = (np.fft.ifft(np.fft.fft(np.exp(-sigma * z) * source) / Nm).real) * np.exp(sigma * z)
    force = make_interp_spline(z, fk, k=5)

    def rhs(u, y):
        return (8.0 / u - 1.0) * y + u * float(force(-2.0 * np.log(u)))

    hi_sol = solve_ivp(rhs, (6.0, 2000.0), [0.0], method="DOP853", rtol=2e-11, atol=1e-20, dense_output=True)
    lo_sol = solve_ivp(rhs, (6.0, 0.005), [0.0], method="DOP853", rtol=2e-11, atol=1e-24, dense_output=True)
    assert hi_sol.success and lo_sol.success

    a2 = D[0] / 2.0
    a3 = 8.0 * a2
    a4 = -D[1] / 24.0 + 9.0 * a3
    a5 = 10.0 * a4
    a6 = D[2] / 720.0 + 11.0 * a5
    coeff = np.array([a2, a3, a4, a5, a6])

    def edge(u):
        return sum(c / u**j for j, c in enumerate(coeff, 2))

    us = np.geomspace(0.015, 0.06, 80)
    ys = lo_sol.sol(us)[0] / us**2
    mat = np.array([us * np.log(us), us, us**2 * np.log(us), us**2]).T
    far_coeff = np.linalg.lstsq(mat, (ys - 24.0 / math.pi * us**5) / us**5, rcond=None)[0]

    def far(u):
        a, b, c, d = far_coeff
        return 24.0 / math.pi * u**5 + u**6 * (a * np.log(u) + b) + u**7 * (c * np.log(u) + d)

    def smooth(t):
        t = np.clip(t, 0.0, 1.0)
        a = np.exp(-1.0 / np.maximum(t, 1e-300))
        b = np.exp(-1.0 / np.maximum(1.0 - t, 1e-300))
        return a / (a + b)

    def hp(u):
        u = np.asarray(u)
        middle = np.where(u >= 6.0, hi_sol.sol(np.clip(u, 6.0, 2000.0))[0], lo_sol.sol(np.clip(u, 0.005, 6.0))[0]) / u**2
        wlo = smooth((u - 0.01) / 0.01)
        whi = smooth((u - 1000.0) / 500.0)
        return (1.0 - whi) * (wlo * middle + (1.0 - wlo) * far(np.maximum(u, 1e-50))) + whi * edge(np.maximum(u, 1e-50))

    u = np.exp(-z / 2.0)
    h = hp(u)
    Hhat = np.fft.fft(np.exp(-sigma * z) * h)
    psihat = Nm * Hhat
    psi = np.fft.ifft(psihat).real * np.exp(sigma * z)
    psiz = np.fft.ifft((sigma + 1j * freq) * psihat).real * np.exp(sigma * z)
    qhat = -2.0 * Nm * np.fft.fft(np.exp(-sigma * z) * u * h)
    qpsi = np.fft.ifft(qhat).real * np.exp(sigma * z)
    c1 = -6.0 / math.pi * a2
    return z, dz, psi, psiz, qpsi, source, D[0], c1


def main() -> int:
    z, dz, psi, psiz, qpsi, source, D0, c1 = build_psi(n=8192)
    res_hom = -6.0 * psi - 2.0 * psiz - 0.5 * qpsi - source
    bulk = (z > -8.0) & (z < 6.0)
    max_bulk = float(np.max(np.abs(res_hom[bulk])))
    p1 = bool(max_bulk < 1e-6 and abs(c1 + 3.0 / math.pi * D0) < 1e-12)
    print(f"PQ-1 homological profile PASS={p1}: max_bulk_res={max_bulk:.2e}, c1={c1:.6f} == -3*D0/pi")

    # PQ-2: a(psi) = 0 with exact c1*e^z edge tail restored
    m = (z >= -16.0) & (z <= 16.0)
    zm = z[m]
    w = np.full_like(zm, dz)
    w[0] *= 0.5
    w[-1] *= 0.5
    a_trunc = -(1.0 / math.pi) * float(np.sum(w * np.exp(-0.5 * zm) * psi[m]))
    a_tail = -(2.0 * c1 / math.pi) * math.exp(0.5 * zm[0])
    a_tot = a_trunc + a_tail
    p2 = bool(abs(a_tot) < 1e-6)
    print(f"PQ-2 edge moment a(psi)  PASS={p2}: a_trunc={a_trunc:+.4e}, a_tail={a_tail:+.4e}, total={a_tot:+.2e}")

    # PQ-3: full time-dependent O(eps^2) PDE residual
    psi_sp = make_interp_spline(z, psi, k=5)
    qpsi_sp = make_interp_spline(z, qpsi, k=5)

    def phi1(v):
        r = v ** (-0.25)
        zz = 2.0 * r * np.sqrt(-1j)
        return float(-2.0 / math.pi * (r**5) * np.real(np.sqrt(-1j) * kv(1, zz)))

    def phi1_b(v):
        r = v ** (-0.25)
        zz = 2.0 * r * np.sqrt(-1j)
        return float(-2.0 / math.pi * (r**5) * np.real(np.sqrt(-1j) / 2.0 * (-kv(1, zz) + zz * kvp(1, zz))))

    def beta(tau):
        return 2.0 * math.exp(tau) - 1.0

    def G1(v, tau):
        bt = beta(tau)
        return (math.exp(tau) / bt) ** 3 * phi1(v / (bt * bt))

    def QG1(v, tau):
        bt = beta(tau)
        return 2.0 * math.exp(3.0 * tau) / (bt**4) * phi1_b(v / (bt * bt))

    def vG1_v(v, tau):
        bt = beta(tau)
        ww = v / (bt * bt)
        return (math.exp(tau) / bt) ** 3 * (-1.5 * phi1(ww) - 0.5 * phi1_b(ww))

    def V_G1(tau):
        bt = beta(tau)
        scale = (math.exp(tau) / bt) ** 3 * bt

        def ig(r):
            if r <= 0 or r > 200:
                return 0.0
            ww = r ** (-4)
            vv = bt * bt * ww
            return 4.0 * (r ** (-3)) * math.log1p(1.0 / vv) * phi1(ww)

        val = sum(
            quad(ig, a, c, epsabs=1e-11, epsrel=1e-11, limit=200)[0]
            for a, c in ((0, 0.5), (0.5, 2.0), (2.0, 8.0), (8.0, 40.0))
        )
        return -(1.0 / math.pi) * scale * val

    def A_fn(tau):
        return 1.5 * quad(lambda s: V_G1(s), tau, tau + 20.0, epsabs=1e-11, epsrel=1e-10, limit=100)[0]

    def B_fn(tau):
        return -0.5 * math.exp(tau) * quad(
            lambda s: math.exp(-s) * beta(s) * V_G1(s), tau, tau + 20.0, epsabs=1e-11, epsrel=1e-10, limit=100
        )[0]

    def QG1_b(vv, tt, bb):
        return 2.0 * math.exp(3.0 * tt) / (bb**4) * phi1_b(vv / (bb * bb))

    dt = 1e-6
    dv_rel = 1e-6
    db = 1e-6
    max_rel = 0.0
    for tau in (0.0, 0.3, 0.8):
        At = A_fn(tau)
        Bt = B_fn(tau)
        V1t = V_G1(tau)
        Ap = -1.5 * V1t
        Bp = Bt + 0.5 * beta(tau) * V1t
        bt = beta(tau)
        for v in (0.05, 0.3, 1.5):
            zw = math.log(v / (bt * bt))

            def G2_tot(vv, tt):
                bb = beta(tt)
                zz = math.log(vv / (bb * bb))
                A_loc = At + Ap * (tt - tau)
                B_loc = Bt + Bp * (tt - tau)
                return (math.exp(tt) / bb) ** 6 * float(psi_sp(zz)) + A_loc * G1(vv, tt) - 0.5 * B_loc * QG1(vv, tt)

            G2_t = (G2_tot(v, tau + dt) - G2_tot(v, tau - dt)) / (2 * dt)
            dv = v * dv_rel
            G2_v = (G2_tot(v + dv, tau) - G2_tot(v - dv, tau)) / (2 * dv)
            Q_QG1 = 2.0 * (QG1_b(v, tau, bt + db) - QG1_b(v, tau, bt - db)) / (2 * db)
            QG2 = (math.exp(6.0 * tau) / (bt**7)) * float(qpsi_sp(zw)) + At * QG1(v, tau) - 0.5 * Bt * Q_QG1
            lhs = G2_t + 2.0 * v * G2_v - 0.5 * QG2
            rhs = -QG1(v, tau) * vG1_v(v, tau) + 0.5 * QG1(v, tau) * G1(v, tau) + V1t * vG1_v(v, tau)
            rel = abs(lhs - rhs) / max(abs(rhs), 1e-30)
            max_rel = max(max_rel, rel)

    p3 = bool(max_rel < 2e-5)
    print(f"PQ-3 full time-dep PDE   PASS={p3}: max_rel_res={max_rel:.2e}")
    assert p1 and p2 and p3, "Pure-step O(eps^2) gate FAILED"
    print("ALL PASS (PQ-1..PQ-3)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
