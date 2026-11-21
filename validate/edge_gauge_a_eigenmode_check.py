#!/usr/bin/env python3
"""Probe: explicit eigenfunctions of the UNCONSTRAINED Gauge-(a) linearization.

Claim under test (reviewer finding, 2026-11-21; see
docs/findings_session_2026_11_21_gauge_a_eigenmodes.md):

  For every lam with Re(lam) > -2*sigma, the function

      G_lam(z) := (1 / 2 pi i) * int_{Re s = s0} e^{s z} X(s; lam) ds,
                  max(0, -Re(lam)/2) < s0 < 1,
      X(s; lam) = Gamma(2s + lam) Gamma(2s + 2) Gamma(s) Gamma(1 - s),

  solves   lam * G = -2 dG/dz + (1/2) Qtilde_op G   on all of R, where
  Qtilde_op G (1+v) = (1/v) int_0^v ( H[g](1+u) - H[g](1) ) du,  g = v^{-1/2} G.

  The lower limit on s0 keeps the contour to the right of the pole of Gamma(2s + lam)
  at s = -lam/2 (s0 = sigma is admissible).  A contour to the left of that pole gives a
  different function (it differs by the residue term A_lam e^{-lam z/2}).

  G_lam is in L^2_sigma(R) whenever max(0, -Re(lam)/2) < sigma < 1 (the upper bound is the
  edge: G ~ r_1 v gives int_0 v^{1 - 2 sigma} dv), is smooth up to the edge with G = O(v),
  and has H[g](1) = -2 Gamma(1 + lam) != 0.  Far field (z -> +inf):
    Re(lam) > 0:                 G_lam -> Gamma(lam);
    -2 sigma < Re(lam) <= 0,     G_lam = A_lam e^{-lam z/2} + Gamma(lam) + o(1),
    lam != 0:                    A_lam = (1/2) Gamma(2-lam) Gamma(-lam/2) Gamma(1+lam/2).
  So the operator -2 d/dz + (1/2) Qtilde_op WITHOUT the domain condition Ghat(1/2) = 0
  has every such lam (including Re lam > 0) as an eigenvalue.  The paper's Theorem 1
  (kernel trivial on Re lam > -2 sigma) holds only on the domain WITH that condition.

  Correction 2026-11-21 (Codex finding, independently verified by Claude): the first
  version of this text allowed any 0 < s0 < 1 and claimed G_lam -> Gamma(lam) for all
  Re(lam) > -2 sigma.  Both are wrong for Re(lam) <= 0; EM-5 below tests lam = -1/4.

What this script checks (exit 0 iff EM-1..EM-5 all PASS; PASS = "finding reproduced"):

  EM-1  Real-space eigen-identity, by quadrature of the Hilbert transform, with NO Mellin
        identity used on the right-hand side:
            lam*G + (lam+4) v G' + 2 v^2 G''  =  (1/2) D(v),
            D(v) = H[g](1+v) - H[g](1) = (v/pi) PV int_0^inf w^{-3/2} G(w) / (v - w) dw.
        (This is d/dv of  v * [lam*G + 2 v G' - (1/2) Qtilde_op G] = 0; both sides vanish
        at v -> 0, so it is equivalent to the eigen-equation.)
  EM-2  Closed forms: H[g](1) = -2 Gamma(1+lam);  G(+inf) = Gamma(lam);  the inversion
        contour is immaterial (Re s = 0.9 vs 0.1).
  EM-3  Edge germ: G ~ r_1 v + r_2 v^2 with r_m = (-1)^{m+1} Gamma(2m+lam) Gamma(2m+2).
        The series has zero radius of convergence, so G_lam is NOT in H^1(D(1,1)):
        no contradiction with the edge-analytic nonlinear theorem.
  EM-4  Discrepancy with the repo's discretized operator (src/stability_edge_z.py, the
        matrices used by edge_nonlinear_comoving_closure.py): the discrete Gauge-(a)
        matrix has no unstable eigenvalue, and marching G_1 with it makes H[g](1) decay
        like e^{-T}, where the continuum mode grows like e^{+T}.  The grid stops at
        z_min, so the edge term is absent from the discretization.
  EM-5  Negative lam (lam = -1/4, pole at s = 1/8): admissible contours agree; a contour
        left of the pole differs by exactly A_lam e^{z/8}; the far field is
        A_lam e^{z/8} + Gamma(lam) + o(1), NOT the constant Gamma(lam); the real-space
        eigen-identity and H[g](1) = -2 Gamma(1+lam) still hold.

EM-1..EM-3 use lam = 0.3, 1, 3 (positive only); EM-5 is the only negative-lam test, and
only real lam is tested.

Not claimed: anything about the nonlinear PDE, the frozen-exponent operator on its
constrained domain, or Conjecture 2.4.  Review form: EM-1..EM-4 single reviewer (Claude);
the contour/far-field correction and EM-5 are an independent Claude verification of
Codex findings.  Nothing here machine-certifies the PDE proof.
"""
from __future__ import annotations

import math
import os
import sys
import warnings

import numpy as np
from scipy.integrate import IntegrationWarning, quad
from scipy.interpolate import CubicSpline
from scipy.special import gamma, loggamma

warnings.simplefilter("ignore", IntegrationWarning)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FAILURES: list[str] = []

X_GRID = np.linspace(-40.0, 40.0, 8001)      # x = ln v
TAU = np.linspace(-12.0, 12.0, 4801)          # |X| ~ e^{-3 pi |tau|}: negligible beyond


def report(name: str, ok: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}  {detail}")
    if not ok:
        FAILURES.append(name)


def X_mult(s, lam):
    return np.exp(loggamma(2 * s + lam) + loggamma(2 * s + 2) + loggamma(s) + loggamma(1 - s))


def inv_mellin(x, lam, s0):
    """Rows: G, v G', v^2 G''  at v = e^x, from (1/2pi) int P(s) e^{s x} X(s) dtau, P = 1, s, s(s-1)."""
    s = s0 + 1j * TAU
    Xs = X_mult(s, lam)
    stack = np.stack([Xs, s * Xs, s * (s - 1) * Xs], axis=1)
    out = np.empty((3, len(x)))
    dt = TAU[1] - TAU[0]
    for i0 in range(0, len(x), 400):
        xb = x[i0:i0 + 400]
        out[:, i0:i0 + 400] = (np.exp(np.outer(xb, s)) @ stack).real.T * dt / (2 * np.pi)
    return out


def far_contour(lam):
    """Far-field contour for real lam: small for conditioning, but right of the pole at s = -lam/2."""
    return max(0.1, -0.5 * lam + 0.15)


def g_lam(x, lam):
    """Contour near Re s = 1 for the edge side (x < 0), near the left end of the admissible
    strip max(0, -lam/2) < Re s < 1 for the far field: best conditioning."""
    out = np.empty((3, len(x)))
    m = x < 0
    if m.any():
        out[:, m] = inv_mellin(x[m], lam, 0.9)
    if (~m).any():
        out[:, ~m] = inv_mellin(x[~m], lam, far_contour(lam))
    return out


def half_D(v, Gs):
    """(1/2) D(v) = (v / 2 pi) PV int_0^inf w^{-3/2} G(w) / (v - w) dw, in real space."""
    f = lambda xx: math.exp(-0.5 * xx) * float(Gs(xx))          # w^{-3/2} G dw = e^{-x/2} G dx
    h = lambda w: -(w ** -1.5) * float(Gs(math.log(w)))
    mid = quad(h, v / 2, 2 * v, weight="cauchy", wvar=v, limit=400)[0]
    lo = quad(lambda xx: f(xx) / (v - math.exp(xx)), X_GRID[0], math.log(v / 2), limit=800)[0]
    hi = quad(lambda xx: f(xx) / (v - math.exp(xx)), math.log(2 * v), X_GRID[-1], limit=800)[0]
    return 0.5 * (v / math.pi) * (lo + mid + hi)


def main() -> int:
    print("Unconstrained Gauge-(a) eigenmode probe (EM-1..EM-5)")
    lams = (0.3, 1.0, 3.0)
    worst_rel, detail1 = 0.0, []
    ok2, detail2 = True, []
    ok3, detail3 = True, []

    for lam in lams:
        G, vG1, v2G2 = g_lam(X_GRID, lam)
        Gs, G1s, G2s = CubicSpline(X_GRID, G), CubicSpline(X_GRID, vG1), CubicSpline(X_GRID, v2G2)

        # ---- EM-1: real-space eigen-identity ----
        rel = 0.0
        for v in (0.05, 0.3, 1.0, 4.0):
            xv = math.log(v)
            lhs = lam * float(Gs(xv)) + (lam + 4) * float(G1s(xv)) + 2 * float(G2s(xv))
            rel = max(rel, abs(lhs - half_D(v, Gs)) / abs(lhs))
        worst_rel = max(worst_rel, rel)
        detail1.append(f"lam={lam}: {rel:.1e}")

        # ---- EM-2: closed forms ----
        a_num = -(1 / math.pi) * quad(lambda xx: math.exp(-0.5 * xx) * float(Gs(xx)), X_GRID[0], X_GRID[-1], limit=800)[0]
        a_exact = -2 * gamma(1 + lam)
        c_hi = inv_mellin(np.array([0.0]), lam, 0.9)[0, 0]
        c_lo = inv_mellin(np.array([0.0]), lam, 0.1)[0, 0]
        ok_here = abs(a_num - a_exact) < 1e-6 * abs(a_exact) and abs(c_hi - c_lo) < 1e-9
        if lam >= 1.0:   # for small lam the e^{-lam z/2} correction is still visible at z = 40
            ok_here = ok_here and abs(G[-1] - gamma(lam)) < 1e-6
        ok2 = ok2 and ok_here
        detail2.append(f"lam={lam}: H[g](1)={a_num:.6f} (-2Gamma(1+lam)={a_exact:.6f}), G(+inf)={G[-1]:.5f} (Gamma(lam)={gamma(lam):.5f})")

        # ---- EM-3: edge germ r_1 v + r_2 v^2, factorially divergent ----
        r1 = gamma(2 + lam) * gamma(4)
        r2 = -gamma(4 + lam) * gamma(6)
        i = int(np.argmin(np.abs(X_GRID + 14.0)))
        v = math.exp(X_GRID[i])
        r1_num = G[i] / v
        r2_num = (G[i] / v - r1) / v
        ok3 = ok3 and abs(r1_num - r1) < 1e-3 * r1 and abs(r2_num - r2) < 2e-2 * abs(r2)
        ratios = [gamma(2 * m + 2 + lam) * gamma(2 * m + 4) / (gamma(2 * m + lam) * gamma(2 * m + 2)) for m in (1, 5, 10)]
        detail3.append(f"lam={lam}: r1={r1_num:.4f} ({r1:.4f}), r2={r2_num:.1f} ({r2:.1f}), |r_(m+1)/r_m| at m=1,5,10: {ratios[0]:.0f}, {ratios[1]:.0f}, {ratios[2]:.0f}")

    report("EM-1 G_lam solves lam*G = -2 G_z + (1/2) Qtilde_op G in real space (rel. residual < 1e-6)",
           worst_rel < 1e-6, "; ".join(detail1))
    report("EM-2 H[g](1) = -2 Gamma(1+lam) != 0, G(+inf) = Gamma(lam), contour-independent",
           ok2, " | ".join(detail2))
    report("EM-3 edge germ r_m = (-1)^(m+1) Gamma(2m+lam) Gamma(2m+2): smooth, O(v), NOT analytic in a disk",
           ok3, " | ".join(detail3))

    # ---- EM-4: the repo's discretized Gauge-(a) operator does not carry these modes ----
    sys.path.insert(0, os.path.join(ROOT, "src"))
    try:
        from scipy.linalg import expm
        from stability_edge_z import build_blocks_z
    except Exception as exc:  # pragma: no cover
        report("EM-4 discrepancy with the repo's discretized operator", False, f"import failed: {exc}")
    else:
        zmin, zmax = -10.0, 10.0
        B = build_blocks_z(N=1000, zmin=zmin, zmax=zmax, beta=0.0, edge_bc=True)
        z, dz, Y, Q_op, Dz = B["z"], B["dz"], B["Y"], B["Q_op"], B["Dz"]
        n = len(z)
        w = np.full(n, dz)
        w[0] *= 0.5
        w[-1] *= 0.5
        h1 = -(1.0 / np.pi) * np.exp(-0.5 * z) * w                       # discrete H[g](1), cut at z_min
        Q_exact = Q_op + np.outer(math.exp(zmin) / (Y - 1.0), h1)        # as in edge_nonlinear_comoving_closure.py
        L = -2.0 * Dz + 0.5 * (Q_exact - np.outer(np.ones(n), h1))
        max_re = float(np.linalg.eigvals(L).real.max())
        G1 = g_lam(z, 1.0)[0]
        a0 = float(h1 @ G1)
        ratios = {T: float(h1 @ (expm(T * L) @ G1)) / a0 for T in (0.5, 1.0, 2.0)}
        decays = all(abs(ratios[T] - math.exp(-T)) < 0.2 * math.exp(-T) for T in ratios)
        report("EM-4 repo's discrete Gauge-(a) matrix: no unstable eigenvalue, and H[g](1) of G_1 decays ~e^(-T) (continuum: e^(+T))",
               max_re < 0 and decays,
               f"max Re(eig)={max_re:.4f}; H[g](1)(T)/H[g](1)(0) at T=0.5,1,2: "
               + ", ".join(f"{ratios[T]:.3f}" for T in ratios)
               + f" vs e^-T: {math.exp(-0.5):.3f}, {math.exp(-1):.3f}, {math.exp(-2):.3f}; continuum e^+T: {math.exp(0.5):.2f}, {math.exp(1):.2f}, {math.exp(2):.2f}")

    # ---- EM-5: negative lam: admissible contour and far-field asymptotics ----
    lam = -0.25
    s_pole = -0.5 * lam                                                   # 1/8
    A = 0.5 * gamma(2 - lam) * gamma(-0.5 * lam) * gamma(1 + 0.5 * lam)   # residue of X at s = -lam/2
    zs = np.array([-6.0, 0.0, 4.0, 16.0])
    g_a = inv_mellin(zs, lam, 0.25)[0]
    g_b = inv_mellin(zs, lam, 0.9)[0]
    g_left = inv_mellin(zs, lam, 0.1)[0]                                  # 0.1 < 1/8: left of the pole
    contours_agree = float(np.max(np.abs(g_a - g_b) / np.abs(g_a))) < 1e-8
    jump_is_residue = float(np.max(np.abs((g_a - g_left) - A * np.exp(s_pole * zs)) / (A * np.exp(s_pole * zs)))) < 1e-8
    zf = np.array([16.0, 24.0, 32.0])
    gf = inv_mellin(zf, lam, 0.25)[0]
    rem = gf - A * np.exp(s_pole * zf) - gamma(lam)                       # should be o(1)
    s_next = -0.5 * (1 + lam)                                             # next pole of Gamma(2s + lam)
    c_next = -0.5 * gamma(2 * s_next + 2) * gamma(s_next) * gamma(1 - s_next)
    far_ok = (abs(rem[0]) > abs(rem[1]) > abs(rem[2])
              and abs(rem[2]) < 1e-4
              and abs(rem[1] - c_next * math.exp(s_next * zf[1])) < 1e-2 * abs(rem[1])
              and abs(gf[2] - gamma(lam)) > 100.0)                        # the constant-limit claim fails
    Gn, vG1n, v2G2n = g_lam(X_GRID, lam)
    Gs, G1s, G2s = CubicSpline(X_GRID, Gn), CubicSpline(X_GRID, vG1n), CubicSpline(X_GRID, v2G2n)
    rel5 = 0.0
    for v in (0.05, 0.3, 1.0, 4.0, 30.0):
        xv = math.log(v)
        lhs = lam * float(Gs(xv)) + (lam + 4) * float(G1s(xv)) + 2 * float(G2s(xv))
        rel5 = max(rel5, abs(lhs - half_D(v, Gs)) / abs(lhs))
    x_max = X_GRID[-1]
    a_tail = -(1 / math.pi) * A * math.exp(-0.5 * (1 + lam) * x_max) * 2 / (1 + lam)   # leading tail beyond the grid
    a_num = -(1 / math.pi) * quad(lambda xx: math.exp(-0.5 * xx) * float(Gs(xx)), X_GRID[0], x_max, limit=800)[0] + a_tail
    a_exact = -2 * gamma(1 + lam)
    a_ok = abs(a_num - a_exact) < 1e-6 * abs(a_exact)
    report("EM-5 lam=-1/4: contour must satisfy s0 > -lam/2; far field is A e^(-lam z/2) + Gamma(lam) + o(1), not Gamma(lam)",
           contours_agree and jump_is_residue and far_ok and rel5 < 1e-6 and a_ok,
           f"A={A:.10f}; G(16)={gf[0]:.10f}, G(32)={gf[2]:.6f} vs Gamma(lam)={gamma(lam):.6f}; "
           f"remainder at z=16,24,32: {rem[0]:.2e}, {rem[1]:.2e}, {rem[2]:.2e}; "
           f"eigen-identity rel. residual {rel5:.1e}; H[g](1)={a_num:.7f} (-2Gamma(1+lam)={a_exact:.7f})")

    if FAILURES:
        print(f"PROBE FAIL: {FAILURES}")
        return 1
    print("PROBE PASS: unconstrained Gauge-(a) eigenmodes reproduced (EM-1..EM-5)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
