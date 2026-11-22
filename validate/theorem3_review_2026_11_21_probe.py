#!/usr/bin/env python3
"""Probe: reproduces the numerically checkable findings of the review of 2026-11-21 of commit
ef0eaea ("Prove coupled (G,h) nonlinear stability with zero derivative loss (Theorem 3, Lean 4,
CH-1..4)"); see docs/findings_session_2026_11_21_theorem3_review.md.

PASS here means "the review finding reproduces", NOT "the manuscript's claim holds".

  Q1  Step 1 of Theorem 3.  The smoothing bound ||d_z^2 (Q_op G)||_{L^2(dz)} <= ||G_z||_{L^2_{1/2}}
      holds (sharp constant 1).  With the manuscript's own definition p = V_G(1) - Q_op G, the
      stated ||p_zz||_{L^2(dz)} <= (1/2) 2^{1/4} sqrt(E(G)) is false; 2^{1/4} sqrt(E(G)) is true.
      Checked by Plancherel on Re w = 1/2 with the exact multiplier for smooth bumps supported
      in z > 0 (so G vanishes on the disk and E = 20 X_far + Y_far).
  Q2  The repo's discrete Q_op (stability_edge_z.build_blocks_z, reused by
      coevolve_moving_frame.build) tends to H[g](1) at the edge.  It is Q_op, not the subtracted
      Qtilde_op = Q_op - H[g](1) that the gate edge_coupled_h_nonlinear_closure.py assumes.
  Q3  CH-4 of that gate, rerun on its own grid, data and time step with two right-hand sides:
      (a) as written in the gate; (b) the Gauge-(a) system of the manuscript (subtract H[g](1),
      p = V_G(1) - Q_op G).  (a) decays; (b) is non-finite before tau = 3.
  Q4  Mechanism of (b).  The discretised Gauge-(a) linear operator
          L_a = M_GG - (1/2) 1 (x) h1        (h1 @ G = H[g](1))
      has one eigenvalue with large positive real part that grows with the edge depth, while
      M_GG has none; and the explicit eigenfunctions G_lambda of OBS-09 satisfy
      L_a G_lambda ~ lambda G_lambda on the grid.  A finite-difference grid cannot impose the
      edge analyticity that Theorem 2 / Theorem 3 use to exclude these modes.
  Q5  The weighted transient of M_GG quoted in the manuscript (peak about 1.17 at gamma = 0.5,
      12.7 at gamma = 0.75, z_min = -14), recomputed with scipy's expm_multiply (no eigen-
      decomposition); and its driver: at the edge G grows initially like (1/2) H[g](1) tau.

Run:  python3 validate/theorem3_review_2026_11_21_probe.py        (about 15 s)
Exit 0 iff Q1..Q5 all reproduce.  The probe does not import the reviewed gate, so it runs with
or without commit ef0eaea.

Review form: one Claude session, independent of the author of ef0eaea; no separate reviewers.
Nothing in this file certifies or refutes Theorem 3.
"""
from __future__ import annotations

import math
import os
import sys
import warnings

import numpy as np
from scipy.sparse.linalg import expm_multiply

warnings.simplefilter("ignore", RuntimeWarning)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "validate"))

from coevolve_moving_frame import build  # noqa: E402
from stability_edge_z import build_blocks_z  # noqa: E402
import edge_gauge_a_eigenmode_check as EM  # noqa: E402

FAILURES: list[str] = []


def report(name: str, ok: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n         {detail}")
    if not ok:
        FAILURES.append(name)


def c_sym(w):
    """Mellin symbol of T_c:  (1/2) Q_op G  <->  c(w) Ghat(w)."""
    return -np.tan(np.pi * w) / (2 * w + 1)


def q1() -> None:
    z = np.linspace(-1.0, 4.0, 250001)
    dz = z[1] - z[0]
    taus = np.linspace(-300.0, 300.0, 12001)
    dt = taus[1] - taus[0]
    s = 1j * taus
    m2_half = s ** 2 * c_sym(s + 0.5) / (s + 0.5)      # d_z^2((1/2) Q_op G)  <->  m2_half * (G_z)^
    m2_half[~np.isfinite(m2_half)] = 0.0               # s = 0: double zero times a simple pole
    rows, planch, r_half, r_E = [], [], [], []
    for (a, b) in ((0.02, 0.22), (0.02, 0.12), (0.5, 1.5)):
        x = (z - a) / (b - a)
        inside = (x > 0) & (x < 1)
        G = np.zeros_like(z)
        G[inside] = np.exp(-1.0 / (x[inside] * (1 - x[inside])))
        Gz = np.gradient(G, dz)
        E = 20 * np.sum(np.exp(-0.5 * z) * G ** 2) * dz + np.sum(np.exp(-0.5 * z) * Gz ** 2) * dz
        n_half = math.sqrt(np.sum(np.exp(-z) * Gz ** 2) * dz)          # ||G_z||_{L^2_{1/2}}
        msk = Gz != 0
        zz, gg = z[msk], Gz[msk]
        F = np.array([np.sum(np.exp(-(0.5 + 1j * t) * zz) * gg) * dz for t in taus])
        planch.append(math.sqrt(np.sum(np.abs(F) ** 2) * dt / (2 * math.pi)) / n_half)
        nQ = math.sqrt(np.sum(np.abs(2 * m2_half * F) ** 2) * dt / (2 * math.pi))   # ||d_z^2 (Q_op G)||
        r_half.append(nQ / n_half)
        r_E.append(nQ / math.sqrt(E))
        rows.append(f"[{a},{b}]: {nQ / n_half:.4f} / {nQ / math.sqrt(E):.4f}")
    stated, true_c = 0.5 * 2 ** 0.25, 2 ** 0.25
    ok = (all(abs(p - 1) < 0.01 for p in planch) and max(r_half) <= 1.0005
          and max(r_E) > stated and max(r_E) <= true_c
          and abs(float(np.max(np.abs(2 * m2_half))) - 1.0) < 1e-3)
    report("Q1 smoothing bound holds with constant 1; the stated (1/2) 2^(1/4) for p_zz does not",
           ok, "||d2 Q_op G|| / ||G_z||_(1/2)  and  / sqrt(E), by support of the bump: " + "; ".join(rows)
           + f".  Stated constant {stated:.4f}, true constant {true_c:.4f}, sup|2 m_2| = "
           f"{float(np.max(np.abs(2 * m2_half))):.5f}")


def q2():
    B = build(N=800, zmin=-10.0, zmax=10.0)
    z, dz, Q = B["z"], B["dz"], B["Q_op"]
    wt = np.full_like(z, dz)
    wt[0] *= 0.5
    wt[-1] *= 0.5
    h1 = -(1.0 / np.pi) * np.exp(-0.5 * z) * wt
    G = 0.15 * np.exp(-((z - 1.0)) ** 2)
    H, QG, V = float(h1 @ G), Q @ G, float(B["b1"] @ G)
    i4, i8 = int(np.argmin(np.abs(z + 4.0))), int(np.argmin(np.abs(z + 8.0)))
    ok = abs(QG[i4] - H) < 0.02 * abs(H)
    report("Q2 the discrete Q_op tends to H[g](1) at the edge (it is not Qtilde_op)", ok,
           f"(Q_op G)(z=-4) = {QG[i4]:+.5f}, H[g](1) = {H:+.5f}; velocity perturbation at z=-8: "
           f"gate's V - Q - H = {V - QG[i8] - H:+.5f}, manuscript's p = V - Q_op G = {V - QG[i8]:+.5f}")
    return B, h1


def q3(B, h1) -> None:
    z, dz, Q_op, Dz, b1, M_GG = B["z"], B["dz"], B["Q_op"], B["Dz"], B["b1"], B["M_GG"]
    n = len(z)
    D_up = np.zeros((n, n))
    for i in range(1, n):
        D_up[i, i], D_up[i, i - 1] = 1.0 / dz, -1.0 / dz

    def rhs(G, q, Bamp, variant):
        V, H = float(b1 @ G), float(h1 @ G)
        QG = Q_op @ G
        Gz, qz = D_up @ G, D_up @ q
        if variant == "gate":                       # as written in edge_coupled_h_nonlinear_closure.py
            Qt, p = QG, V - QG - H
            dG = M_GG @ G + Bamp * q + p * Gz + 0.5 * Qt * G
        else:                                       # manuscript: Qtilde = Q_op - H, p = V - Q_op G
            Qt, p = QG - H, V - QG
            dG = M_GG @ G - 0.5 * H + Bamp * q + p * Gz + 0.5 * Qt * G
        dq = (-2.0 + p) * qz - q + (0.5 * p - (Dz @ QG)) * q
        dG[0] = 0.0
        dq[0] = 0.0
        return dG, dq, -Bamp * H

    def e_tot(G, q):
        m = z >= -math.log(2.0)
        w = np.exp(-0.5 * z[m]) * dz
        return float(np.sum(w * (20 * G[m] ** 2 + (Dz @ G)[m] ** 2)) + 100.0 * np.sum(w * (q[m] ** 2 + (Dz @ q)[m] ** 2)))

    T_end = 3.0
    nsteps = int(round(T_end / (0.2 * dz)))
    dt = T_end / nsteps
    out, ok = [], True
    for variant in ("gate", "gauge_a"):
        for Ga in (0.0, 0.15):
            G = Ga * np.exp(-((z - 1.0)) ** 2)
            q = 0.2 * np.exp(-((z - 1.5)) ** 2)
            G[0] = q[0] = 0.0
            Bamp, E0, t_bad = 1.0, None, None
            E0 = e_tot(G, q)
            for k in range(nsteps):
                k1 = rhs(G, q, Bamp, variant)
                k2 = rhs(G + 0.5 * dt * k1[0], q + 0.5 * dt * k1[1], Bamp + 0.5 * dt * k1[2], variant)
                k3 = rhs(G + 0.5 * dt * k2[0], q + 0.5 * dt * k2[1], Bamp + 0.5 * dt * k2[2], variant)
                k4 = rhs(G + dt * k3[0], q + dt * k3[1], Bamp + dt * k3[2], variant)
                G = G + (dt / 6) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
                q = q + (dt / 6) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
                Bamp = Bamp + (dt / 6) * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2])
                if not np.all(np.isfinite(G)) or float(np.max(np.abs(G))) > 1e6:
                    t_bad = (k + 1) * dt
                    break
            if variant == "gate":
                ratio = e_tot(G, q) / E0
                ok = ok and t_bad is None and ratio < math.exp(-T_end / 16.0)
                out.append(f"gate rhs, G0 amp {Ga}: E_tot ratio {ratio:.2e}")
            else:
                ok = ok and t_bad is not None
                out.append(f"Gauge-(a) rhs, G0 amp {Ga}: " + (f"|G| > 1e6 at tau = {t_bad:.3f}" if t_bad else "stays bounded"))
    report("Q3 CH-4 decays with the gate's right-hand side and blows up with the manuscript's", ok, "; ".join(out))


def q4() -> None:
    tops, det, ok = [], [], True
    for (N, zmin) in ((800, -10.0), (960, -12.0)):
        B = build(N=N, zmin=zmin, zmax=10.0)
        z, dz, M = B["z"], B["dz"], B["M_GG"]
        wt = np.full_like(z, dz)
        wt[0] *= 0.5
        wt[-1] *= 0.5
        h1 = -(1.0 / np.pi) * np.exp(-0.5 * z) * wt
        La = M - 0.5 * np.outer(np.ones_like(z), h1)
        mt, ea = float(np.linalg.eigvals(M).real.max()), np.linalg.eigvals(La)
        tops.append(float(ea.real.max()))
        ok = ok and mt < 0 and tops[-1] > 10.0
        det.append(f"z_min={zmin:+.0f}: max Re eig M_GG = {mt:+.4f}, max Re eig L_a = {tops[-1]:+.2f} "
                   f"({int(np.sum(ea.real > 1e-8))} with Re > 0)")
        if N == 800:
            mid = (z > zmin + 3.0) & (z < 6.0)
            for lam in (0.5, 1.0, 2.0):
                Gl = EM.g_lam(z, lam)[0]
                ra = float(np.max(np.abs((La @ Gl - lam * Gl)[mid])) / np.max(np.abs(lam * Gl[mid])))
                rt = float(np.max(np.abs((M @ Gl - lam * Gl)[mid])) / np.max(np.abs(lam * Gl[mid])))
                ok = ok and ra < 0.1 and rt > 0.5
                det.append(f"lambda={lam}: residual of G_lambda under L_a {ra:.3f}, under M_GG {rt:.2f}")
    ok = ok and tops[1] > 1.5 * tops[0]
    report("Q4 the discretised Gauge-(a) operator has a growing mode that M_GG lacks; the OBS-09 modes fit it", ok, "; ".join(det))


def q5() -> None:
    B = build_blocks_z(N=1400, zmin=-14.0, zmax=12.0, edge_bc=True)
    z, M = B["z"], B["M_GG"]
    dz = z[1] - z[0]
    G0 = np.exp(-((z - 1.0)) ** 2)
    H = -(1.0 / math.pi) * float(np.sum(np.exp(-0.5 * z) * G0) * dz)
    tt = np.linspace(0.0, 0.5, 201)
    Gs = expm_multiply(M, G0, start=0.0, stop=0.5, num=201, endpoint=True)
    peaks = {}
    for g in (0.5, 0.75):
        w = np.exp(-2 * g * z)
        r = np.sqrt(np.array([np.sum(w * Gs[k] ** 2) for k in range(len(tt))]) / np.sum(w * G0 ** 2))
        peaks[g] = (float(r.max()), float(tt[int(np.argmax(r))]))
    i10, k1 = int(np.argmin(np.abs(z + 10.0))), 1                    # tau = 0.0025
    edge, pred = float(Gs[k1][i10]), 0.5 * H * float(tt[k1])
    ok = (11.5 < peaks[0.75][0] < 14.0 and 1.10 < peaks[0.5][0] < 1.25
          and peaks[0.75][1] < 0.05 and peaks[0.5][1] < 0.05 and abs(edge / pred - 1.0) < 0.15)
    report("Q5 weighted transient peaks of M_GG reproduce with an independent propagator; the edge value drives them", ok,
           f"gamma=0.5: peak {peaks[0.5][0]:.3f} at tau={peaks[0.5][1]:.4f}; gamma=0.75: peak {peaks[0.75][0]:.3f} at "
           f"tau={peaks[0.75][1]:.4f}; G(z=-10, tau=0.0025) = {edge:+.3e} against (1/2) H[g](1) tau = {pred:+.3e}")


def main() -> int:
    print("Theorem-3 review probe 2026-11-21 (Q1..Q5; PASS = review finding reproduces)")
    q1()
    B, h1 = q2()
    q3(B, h1)
    q4()
    q5()
    if FAILURES:
        print(f"PROBE FAIL: {FAILURES}")
        return 1
    print("PROBE PASS: review findings Q1..Q5 reproduce")
    return 0


if __name__ == "__main__":
    sys.exit(main())
