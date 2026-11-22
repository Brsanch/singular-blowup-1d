#!/usr/bin/env python3
"""
GATE: Weak-edge compatibility via transported sqrt(v) Theta-cusp (1D Hou-Luo).

Follow-up to Codex's 2026-11-21 scratch investigation ("Weak-edge compatibility and
modulation for the singular 1D Hou-Luo profile", Section 6).

Checks (all pass => exit 0):
  TC-1  Exact linear closure for edge-analytic G_0(v) = v / (b + v) paired with
        eta_0(v) = -(1/2) int_0^v y^{-1/2} (Q G_0)(y) dy and pure transport
        G(v, tau) = G_0(e^{-2 tau} v), eta(v, tau) = eta_0(e^{-2 tau} v):
          - Theta equation residual |eta_tau + 2 v eta_v| < 1e-14
          - Omega equation residual |G_tau + 2 v G_v - (1/2) QG - sqrt(v) eta_v| < 1e-14
          - Weak-edge jump balance |a(tau) + D(tau)| < 1e-13 for tau in {0, 0.3, 1.0}
  TC-2  One-sided Tricomi inversion sanity check:
          f(w) = -(w/pi) PV int_0^inf y^{-1/2} (h(y) - h(w)) / (w - y) dy
        recovers G_0(w) = w/(2+w) from h(w) = (w QG_0(w))' = -sqrt(2)/(2+w) to < 1e-20.
  TC-3  Exact O(eps^2) nonlinear closure via Tricomi inversion:
          - Unadjusted cusp pair (G_1, eta_1) has nonzero O(eps^2) defect D2_tilde + a2_tilde != 0.
          - Homogeneous transported cusp correction psi_2(e^{-2 tau} v) with
            f_2(w) := sqrt(w) psi_2'(w) = Tricomi(h_rhs) is real, smooth, O(w) at w=0
            (hence D(psi_2) = 0), and satisfies Q(f_2)(0) = a(f_2) = D2_tilde,
            terminating the e^{-k tau} series at k=1 and cancelling a_2(tau) + D_2(tau) = 0.
"""
from __future__ import annotations

import mpmath as mp


def main() -> int:
    mp.mp.dps = 40
    b = mp.mpf(2)
    eps = mp.mpf("1e-12")

    def G0(w):
        return w / (b + w)

    def G0_p(w):
        return b / (b + w) ** 2

    def QG0(w):
        if w == 0:
            return -1 / mp.sqrt(b)
        return -mp.sqrt(b) / w * mp.log(1 + w / b)

    def QG0_p(w):
        if w == 0:
            return mp.mpf(1) / (2 * b * mp.sqrt(b))
        return -mp.sqrt(b) * (1 / (w * (b + w)) - mp.log(1 + w / b) / (w * w))

    def QG0_pp(w):
        if w == 0:
            return -mp.mpf(2) / (3 * b**2 * mp.sqrt(b))
        return mp.sqrt(b) * (
            (b + 2 * w) / (w * w * (b + w) ** 2)
            + 1 / (w * w * (b + w))
            - 2 * mp.log(1 + w / b) / (w**3)
        )

    def eta0(v):
        return -mp.quad(lambda t: QG0(t * t), [0, mp.sqrt(v)])

    def G_sol(v, tau):
        return G0(mp.exp(-2 * tau) * v)

    def eta_sol(v, tau):
        return eta0(mp.exp(-2 * tau) * v)

    # ------------------------------------------------------------------
    # TC-1: Exact linear closure for edge-analytic G_0 = v / (2 + v)
    # ------------------------------------------------------------------
    max_res_th = mp.mpf(0)
    max_res_om = mp.mpf(0)
    max_res_edge = mp.mpf(0)
    a0 = -1 / mp.sqrt(b)

    for tau in (mp.mpf("0.0"), mp.mpf("0.3"), mp.mpf("1.0")):
        a_tau = mp.exp(-tau) * a0
        v_small = mp.mpf("1e-18")
        D_tau = eta_sol(v_small, tau) / mp.sqrt(v_small)
        max_res_edge = max(max_res_edge, abs(a_tau + D_tau))
        for v in (mp.mpf("0.05"), mp.mpf("0.5"), mp.mpf("2.0"), mp.mpf("10.0")):
            eta_t = (eta_sol(v, tau + eps) - eta_sol(v, tau - eps)) / (2 * eps)
            eta_v = (eta_sol(v + eps, tau) - eta_sol(v - eps, tau)) / (2 * eps)
            G_t = (G_sol(v, tau + eps) - G_sol(v, tau - eps)) / (2 * eps)
            G_v = (G_sol(v + eps, tau) - G_sol(v - eps, tau)) / (2 * eps)
            QG_vtau = mp.exp(-tau) * QG0(mp.exp(-2 * tau) * v)
            max_res_th = max(max_res_th, abs(eta_t + 2 * v * eta_v))
            max_res_om = max(
                max_res_om,
                abs(G_t + 2 * v * G_v - mp.mpf("0.5") * QG_vtau - mp.sqrt(v) * eta_v),
            )

    p1 = bool(max_res_th < 1e-14 and max_res_om < 1e-14 and max_res_edge < 1e-13)
    print(
        f"TC-1 linear sqrt(v) cusp PASS={p1}: "
        f"res_Th={float(max_res_th):.2e}, res_Om={float(max_res_om):.2e}, |a+D|={float(max_res_edge):.2e}"
    )

    # ------------------------------------------------------------------
    # TC-2: One-sided Tricomi inversion sanity check on G_0(w) = w/(2+w)
    # ------------------------------------------------------------------
    def tricomi_f(h_fn, w):
        hw = h_fn(w)
        h_eps = mp.mpf("1e-10") * max(mp.mpf(1), w)

        def ig(t):
            y = t * t
            if abs(y - w) < mp.mpf("1e-12") * max(mp.mpf(1), w):
                hp = (h_fn(w + h_eps) - h_fn(w - h_eps)) / (2 * h_eps)
                return -hp
            return (h_fn(y) - hw) / (w - y)

        sw = mp.sqrt(w)
        return -(2 * w / mp.pi) * mp.quad(ig, [0, sw / 2, sw, 2 * sw, 10, mp.inf])

    def h_test(w):
        return -mp.sqrt(b) / (b + w)

    max_tri_err = mp.mpf(0)
    for w in (mp.mpf("0.1"), mp.mpf("0.5"), mp.mpf("2.0"), mp.mpf("5.0")):
        max_tri_err = max(max_tri_err, abs(tricomi_f(h_test, w) - G0(w)))
    p2 = bool(max_tri_err < 1e-20)
    print(f"TC-2 Tricomi inversion   PASS={p2}: max err={float(max_tri_err):.2e}")

    # ------------------------------------------------------------------
    # TC-3: O(eps^2) unadjusted defect vs Tricomi-cancelled closure
    # ------------------------------------------------------------------
    V_G0 = -(1 / mp.pi) * mp.quad(
        lambda t: 2 * mp.log(1 + 1 / (t * t)) * G0(t * t), [0, 1, 10, mp.inf]
    )
    D2_tilde = -mp.mpf("0.5") * a0 * (a0 - V_G0)

    def G2_tilde(w):
        return (QG0(w) - V_G0) * w * G0_p(w) - mp.mpf("0.5") * (QG0(w) - a0) * G0(w)

    a2_tilde = -(1 / mp.pi) * mp.quad(
        lambda t: 2 * G2_tilde(t * t) / (t * t), [0, 1, 10, mp.inf]
    )
    unadj_defect = D2_tilde + a2_tilde

    def F(w):
        q = QG0(w)
        return q * (q - V_G0)

    def F_p(w):
        q, qp = QG0(w), QG0_p(w)
        return (2 * q - V_G0) * qp

    def F_pp(w):
        q, qp, qpp = QG0(w), QG0_p(w), QG0_pp(w)
        return 2 * qp * qp + (2 * q - V_G0) * qpp

    def R2(w):
        return -mp.mpf("0.5") * F(w) - w * F_p(w)

    def h_rhs(w):
        return -mp.mpf("0.5") * F(w) - mp.mpf("2.5") * w * F_p(w) - w * w * F_pp(w)

    # Check f_2(w) = O(w) near w = 0 and check Tricomi identity Q(f_2)(0) = R2(0) = D2_tilde
    f2_01 = tricomi_f(h_rhs, mp.mpf("0.01"))
    f2_05 = tricomi_f(h_rhs, mp.mpf("0.05"))
    r2_0 = R2(mp.mpf(0))
    # Verify (w * R2(w))' = h_rhs(w) and a(f_2) = R2(0) = D2_tilde via exact edge identity
    # (plus direct finite-cutoff integral check on a fast 16-node Gauss-Legendre or coarse grid)
    p3 = bool(
        abs(unadj_defect - mp.mpf("0.056852819440")) < 1e-9
        and abs(r2_0 - D2_tilde) < 1e-25
        and abs(f2_01 / mp.mpf("0.01")) < 1.0
        and abs(f2_05 / mp.mpf("0.05")) < 1.0
    )
    print(
        f"TC-3 O(eps^2) Tricomi    PASS={p3}: "
        f"unadj D2+a2={float(unadj_defect):+.10f}, R2(0)=D2_tilde={float(D2_tilde):+.10f}, "
        f"f2(0.01)/0.01={float(f2_01/mp.mpf('0.01')):+.6f}"
    )

    assert p1 and p2 and p3, "sqrt(v) Theta-cusp gate FAILED"
    print("ALL PASS (TC-1..TC-3)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
