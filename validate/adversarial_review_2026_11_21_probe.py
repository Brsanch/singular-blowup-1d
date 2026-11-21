#!/usr/bin/env python3
"""Probe: reproduces the three numerically checkable findings of the adversarial review of
2026-11-21 (docs/findings_session_2026_11_21_adversarial_review.md).

PASS here means "the review finding reproduces", NOT "the manuscript's claim holds".

  R1  Weighted-norm transient of the paper's own linear operator M_GG (inflow BC on).
      The manuscript quotes "peak ||.||_w/||.||_0 = 1.000".  The gate nonlinear_shape_march.py
      computes that peak in the UNWEIGHTED norm (gamma = 0) at seven coarse times.  On a fine
      time grid, in L^2(e^{-2 gamma z} dz): gamma = 0.5 peaks near 1.17 and gamma = 0.75 near
      12.7 on the gate's grid; the peak is N-independent and grows with the edge depth z_min.
  R2  Discrete bulk semigroup of transient_growth_edge_conditioned.py in its W-norm:
      ||exp(M_b tau)||_W exceeds 100 by tau = 6 and is still growing at tau = 8.
  R3  Edge term in d/dtau H[g](1).  For the linear Gauge-(a) flow
          d/dtau H[g](1) = -H[g](1) - c_1/3,      c_1 = G_v(v = 0),
      equivalently  int_0^inf v^{-3/2} Qtilde_op G (1+v) dv = (2 pi / 3) c_1.
      So {H[g](1) = 0} is not invariant.  Checked in real space for three functions (two
      with c_1 = 1, one with c_1 = 0) by quadrature of the Hilbert transform, regularised
      with the profile identity H[(Y-1)^{-1/2}] = 0; no Mellin identity is used.

Run:  python3 validate/adversarial_review_2026_11_21_probe.py      (about one minute)
Exit 0 iff R1..R3 all reproduce.

Review form: R1 and R2 were found by an independent Claude reviewer and re-run by the
coordinating Claude session; R3 was found by Codex, extended by an independent Claude reviewer,
and re-derived here.  Nothing in this file certifies or refutes the PDE theorems.
"""
from __future__ import annotations

import math
import os
import sys
import warnings

import numpy as np
from scipy.integrate import IntegrationWarning, quad
from scipy.linalg import expm

warnings.simplefilter("ignore", IntegrationWarning)
warnings.simplefilter("ignore", RuntimeWarning)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "validate"))

from stability_edge_z import build_blocks_z  # noqa: E402

FAILURES: list[str] = []


def report(name: str, ok: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}  {detail}")
    if not ok:
        FAILURES.append(name)


def weighted_peaks(N: int, zmin: float):
    """sup over a fine tau grid of ||e^{M tau} G0||_gamma / ||G0||_gamma for the gate's bump."""
    B = build_blocks_z(N=N, zmin=zmin, zmax=12.0, edge_bc=True)
    M, z = B["M_GG"], B["z"]
    dz = z[1] - z[0]
    G0 = np.exp(-((z - 1.0) / 1.0) ** 2)
    vals, V = np.linalg.eig(M)
    c = np.linalg.solve(V, G0.astype(complex))
    fine = np.concatenate([np.linspace(0.0, 0.2, 81), np.linspace(0.2, 3.0, 113)])
    coarse = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0])        # the gate's sampling
    out = {}
    for g in (0.0, 0.5, 0.75):
        w = np.exp(-2 * g * z)
        n0 = math.sqrt(float(np.sum(w * G0 ** 2) * dz))
        ratio = lambda t: math.sqrt(float(np.sum(w * np.abs(V @ (np.exp(vals * t) * c)) ** 2) * dz)) / n0
        out[g] = (max(ratio(t) for t in fine), max(1.0, max(ratio(t) for t in coarse)))
    return out


def main() -> int:
    print("Adversarial-review probe 2026-11-21 (R1..R3; PASS = review finding reproduces)")

    # ---- R1: weighted-norm transient peak --------------------------------------------------
    rows = []
    for (N, zmin) in ((1000, -10.0), (1400, -14.0), (2000, -14.0), (1600, -16.0)):
        rows.append((N, zmin, weighted_peaks(N, zmin)))
    p75 = {(N, zm): r[0.75][0] for N, zm, r in rows}
    p50 = {(N, zm): r[0.5][0] for N, zm, r in rows}
    gate_metric = max(r[0.0][1] for _, _, r in rows)          # unweighted, coarse times, floored at 1
    ok1 = (
        abs(gate_metric - 1.0) < 5e-3                         # the quoted 1.000 is this number
        and p75[(1400, -14.0)] > 10.0 and p50[(1400, -14.0)] > 1.1
        and p75[(1000, -10.0)] < p75[(1400, -14.0)] < p75[(1600, -16.0)]   # grows with edge depth
        and abs(p75[(1400, -14.0)] - p75[(2000, -14.0)]) < 0.05 * p75[(1400, -14.0)]  # N-independent
    )
    report(
        "R1 weighted transient peak of M_GG: gamma=0.75 peak >> 1, grows with edge depth, N-independent",
        ok1,
        "gate metric (gamma=0, coarse) = %.3f; fine-grid peak (gamma=0.5 / 0.75): " % gate_metric
        + "; ".join(f"N={N},zmin={zm:+.0f}: {r[0.5][0]:.3f} / {r[0.75][0]:.2f}" for N, zm, r in rows),
    )

    # ---- R2: discrete bulk semigroup of the edge-conditioned gate ---------------------------
    import transient_growth_edge_conditioned as T

    M, Dz, Qop, z, dz = T.build_bounded(1000, -12, 8)
    keep = z > -8.0
    zk = z[keep]
    pp = 2.0 * T.GAMMA * 0.5 * (1.0 + np.tanh((zk - T.Z0) / T.S))
    sw = np.exp(-0.5 * (np.cumsum(pp) * dz))
    Mb = -2.0 * Dz[np.ix_(keep, keep)] + 0.5 * Qop[np.ix_(keep, keep)]
    Mt = (sw[:, None] * Mb) / sw[None, :]
    mu_disc = float(np.linalg.eigvalsh(0.5 * (Mt + Mt.T))[-1])
    norms = {t: float(np.linalg.norm(expm(Mt * t), 2)) for t in (0.001, 2.0, 6.0, 8.0)}
    ok2 = norms[6.0] > 100.0 and norms[8.0] > 1.5 * norms[6.0] and mu_disc > 30.0
    report(
        "R2 ||exp(M_b tau)||_W of the edge-conditioned gate (N=1000, z_min=-12, z_cut=-8) exceeds 100 and keeps growing",
        ok2,
        f"mu_disc={mu_disc:+.2f}; ||exp(M_b tau)||_W at tau=0.001, 2, 6, 8: "
        + ", ".join(f"{norms[t]:.3f}" for t in (0.001, 2.0, 6.0, 8.0)),
    )

    # ---- R3: edge term in d/dtau H[g](1) -----------------------------------------------------
    def edge_integral(phi):
        """For G(1+w) = w * phi(w):
        int_0^inf v^{-3/2} Qtilde_op G (1+v) dv = (2/3) int_0^inf u^{-3/2} D(u) du, where
        D(u) = H[g](1+u) - H[g](1) = (u/pi) PV int_0^inf w^{-1/2} phi(w) / (u - w) dw.
        Because PV int_0^inf w^{-1/2} / (u - w) dw = 0 for u > 0 (the profile identity
        H[(Y-1)^{-1/2}] = 0), D(u) = (u/pi) int_0^inf w^{-1/2} (phi(w) - phi(u)) / (u - w) dw,
        a regular integral.  (The direct principal-value form loses all accuracy at small u,
        where two large pieces cancel and the outer weight u^{-3/2} amplifies the error.)"""
        XLO, XHI = -60.0, 60.0

        def D(u):
            xu, pu = math.log(u), phi(u)
            def integrand(x):                                   # w = e^x, w^{-1/2} dw = e^{x/2} dx
                w = math.exp(x)
                if abs(w - u) < 1e-12 * u:
                    return 0.0
                return math.exp(0.5 * x) * (phi(w) - pu) / (u - w)
            pts = sorted({XLO, XHI, xu, *[p for p in (-10.0, -4.0, -1.5, 0.0, 1.5, 4.0)]})
            return (u / math.pi) * sum(quad(integrand, pts[i], pts[i + 1], limit=200)[0] for i in range(len(pts) - 1))

        # outer integral in x = ln u; the tails beyond |x| = 30 are below 1e-6
        outer = lambda x: math.exp(-0.5 * x) * D(math.exp(x))
        val = sum(quad(outer, a, b, limit=200)[0] for a, b in ((-30.0, -6.0), (-6.0, 0.0), (0.0, 6.0), (6.0, 30.0)))
        return (2.0 / 3.0) * val

    tests = [                                                   # (name of G, phi = G/w, c_1 = phi(0))
        ("v e^{-v}", lambda w: math.exp(-w), 1.0),
        ("v^2 e^{-v}", lambda w: w * math.exp(-w), 0.0),
        ("(v - 2 v^2) e^{-2v}", lambda w: (1.0 - 2.0 * w) * math.exp(-2.0 * w), 1.0),
    ]
    ok3, det = True, []
    for name, G, c1 in tests:
        got = edge_integral(G)
        want = 2 * math.pi / 3 * c1
        ok3 = ok3 and abs(got - want) < 2e-5
        det.append(f"{name}: {got:+.6f} (2pi/3*c1 = {want:+.6f})")
    report("R3 int v^{-3/2} Qtilde_op G dv = (2 pi/3) c_1, so d/dtau H[g](1) = -H[g](1) - c_1/3", ok3, "; ".join(det))

    if FAILURES:
        print(f"PROBE FAIL: {FAILURES}")
        return 1
    print("PROBE PASS: review findings R1..R3 reproduce")
    return 0


if __name__ == "__main__":
    sys.exit(main())
