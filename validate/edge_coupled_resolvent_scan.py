#!/usr/bin/env python3
"""
GATE: Coupled (G,h) resolvent scan — does M_Gh destroy the edge-annihilation gain?

OBS-07 scanned pure-Ω ladder resolvent (continuum): edge forcing gain ∝ e^{z0/2}.
This gate asks the *next* open question on the discrete edge-resolved operator:

  For the block-triangular frozen linearization
      [λ−M_GG, −M_Gh;  0, λ−M_hh] [G;h] = [F; f],
  measure weighted gains for three forcings, as the force is translated edge-ward:

  (A) pure-Ω:   F = bump(z−z0), f = 0     — should reproduce gain∝ε class
  (B) pure-Θ:   F = 0, f = bump(z−z0)     — NEW: h drives G only through M_Gh
  (C) coupled:  F = f = bump(z−z0)

Questions:
  - Is max Re spectrum of full M still ≤ 0 (gap intact under coupling)?
  - Does pure-Θ edge drive produce G-response that *grows* as z0→−∞ (M_Gh ~ e^{−z/2}
    blows at edge), or is it still controlled?
  - Compare gain_B(z0) to gain_A(z0): ratio and log-slopes.

Win: new gain table for coupled/Θ-drive not already in OBS-07 (pure Ω continuum).
Kill if pure-Ω only and matches prior within 5% — this gate always includes (B)(C).

Exit 0 iff structural PASS (finite gains, spectrum, and key numbers printed).
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from stability_edge_z import build_blocks_z, spectrum  # noqa: E402

FAILURES: list[str] = []


def report(name: str, ok: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}  {detail}")
    if not ok:
        FAILURES.append(name)


def wnorm(z, dz, F, sigma):
    return float(np.sqrt(np.sum(np.exp(-2.0 * sigma * z) * np.abs(F) ** 2) * dz))


def bump(z, z0, width=1.0):
    return np.exp(-np.cosh((z - z0) / width) + 1.0)  # peaks ~1 at z0, decays both sides


def solve_resolvent(Mgg, Mhh, Mgh, F, f, lam):
    """
    Block-triangular solve:
      (λ − M_hh) h = f
      (λ − M_GG) G = F + M_Gh h
    """
    n = Mgg.shape[0]
    I = np.eye(n)
    # h
    h = np.linalg.solve(lam * I - Mhh, f)
    # G
    rhs_G = F + Mgh @ h
    G = np.linalg.solve(lam * I - Mgg, rhs_G)
    return G, h


def main() -> int:
    sigma = 0.4
    lam = 0.3  # dist to line −0.8 is 1.1
    z0s = np.array([-2.0, -4.0, -6.0, -8.0, -10.0])
    N = 900
    zmin, zmax = -14.0, 8.0

    print("=" * 72)
    print("  Coupled (G,h) resolvent scan (discrete edge-resolved operator)")
    print(f"  σ={sigma}  λ={lam}  N={N}  z∈[{zmin},{zmax}]")
    print("=" * 72)

    B = build_blocks_z(N=N, zmin=zmin, zmax=zmax, edge_bc=True)
    z, dz = B["z"], B["dz"]
    Mgg, Mhh, Mgh = B["M_GG"], B["M_hh"], B["M_Gh"]
    n = len(z)

    # Spectrum of full coupled operator
    Mfull = np.block([[Mgg, Mgh], [np.zeros((n, n)), Mhh]])
    vfull, _ = spectrum(Mfull)
    max_re = float(vfull.real.max())
    npos = int((vfull.real > 1e-6).sum())
    print(f"\n  full coupled max Re eig = {max_re:+.5f}  (#Re>0 = {npos})")
    report("spectrum: no Re>0 on coupled operator", npos == 0 and max_re < 0.05,
           f"maxRe={max_re:+.5f}")

    rows = []
    print(f"\n  {'z0':>5}  {'gain_Ω':>10}  {'gain_Θ→G':>10}  {'gain_cpl':>10}  "
          f"{'||G||_Θ':>10}  {'||h||_Θ':>10}")
    for z0 in z0s:
        Fb = bump(z, z0)
        # (A) pure Ω
        Ga, ha = solve_resolvent(Mgg, Mhh, Mgh, Fb, np.zeros(n), lam)
        nF = wnorm(z, dz, Fb, sigma)
        gA = wnorm(z, dz, Ga, sigma) / nF
        # (B) pure Θ
        Gb, hb = solve_resolvent(Mgg, Mhh, Mgh, np.zeros(n), Fb, lam)
        gB_G = wnorm(z, dz, Gb, sigma) / nF   # G response / ||f||
        gB_h = wnorm(z, dz, hb, sigma) / nF
        # (C) coupled equal force
        Gc, hc = solve_resolvent(Mgg, Mhh, Mgh, Fb, Fb, lam)
        # combined weighted norm: ||G||_σ² + ||h||_σ² under equal force components
        n_state = np.sqrt(wnorm(z, dz, Gc, sigma) ** 2 + wnorm(z, dz, hc, sigma) ** 2)
        gC = n_state / (np.sqrt(2.0) * nF)
        rows.append((z0, gA, gB_G, gC, gB_G, gB_h))
        print(f"  {z0:+5.1f}  {gA:10.4g}  {gB_G:10.4g}  {gC:10.4g}  "
              f"{gB_G:10.4g}  {gB_h:10.4g}")

    z0a = np.array([r[0] for r in rows])
    gA = np.array([r[1] for r in rows])
    gB = np.array([r[2] for r in rows])
    gC = np.array([r[3] for r in rows])

    # log-slopes d log(gain) / d z0  (OBS-07: → +1/2 for pure Ω edge annihilation)
    def log_slope(z0, g):
        return float(np.polyfit(z0, np.log(np.maximum(g, 1e-30)), 1)[0])

    slope_A = log_slope(z0a, gA)
    slope_B = log_slope(z0a, gB)
    slope_C = log_slope(z0a, gC)
    print(f"\n  log-slopes d log(gain)/dz0:  Ω={slope_A:+.4f}  Θ→G={slope_B:+.4f}  "
          f"cpl={slope_C:+.4f}   (OBS-07 pure-Ω → +0.5)")

    # Finite gains
    report("all gains finite and positive", np.all(np.isfinite(gA)) and np.all(gA > 0)
           and np.all(np.isfinite(gB)) and np.all(gB > 0),
           f"gA={gA} gB={gB}")

    # Pure-Ω should decay toward edge roughly like e^{z0/2}: slope near +0.5
    # (allow loose band — discrete truncation + weight)
    report(
        "pure-Ω edge slope in (0.25, 0.75) — annihilation class",
        0.25 < slope_A < 0.75,
        f"slope_A={slope_A:.4f}",
    )

    # KEY NEW: does Θ→G blow up at edge? M_Gh ~ e^{-z/2} suggests possible growth.
    # Measured: mild growth (slope ~ −0.1) saturating to O(1) — NOT e^{-z0/2} blowup.
    # Ratio (Θ→G)/Ω grows large at depth ONLY because pure-Ω annihilates (denom→0);
    # that is expected and not a coupling instability.
    ratio_deep = float(gB[-1] / gA[-1])
    ratio_bulk = float(gB[0] / gA[0])
    print(f"  ratio (Θ→G)/Ω : bulk z0=-2 → {ratio_bulk:.4g};  deep z0=-10 → {ratio_deep:.4g}")
    print(f"  (large deep ratio is Ω-annihilation, not Θ blowup: gThG_deep={gB[-1]:.4g})")

    report(
        "Θ→G absolute gain bounded (gThG < 3 at all depths; no e^{-z0/2} blowup)",
        float(np.max(gB)) < 3.0,
        f"max gThG={float(np.max(gB)):.4g}  slope_B={slope_B:.4f}",
    )

    # Coupled combined gain also O(1)
    report(
        "coupled combined gain O(1) (max gC < 3)",
        float(np.max(gC)) < 3.0,
        f"gC range [{float(np.min(gC)):.4g}, {float(np.max(gC)):.4g}]",
    )

    # Structural contrast: Ω dies, Θ→G saturates — slopes have opposite sign
    report(
        "Ω annihilates while Θ→G saturates (slope_A>0.25 and gThG_deep/gThG_bulk < 4)",
        slope_A > 0.25 and float(gB[-1] / gB[0]) < 4.0,
        f"slope_A={slope_A:.4f}  gThG growth factor={float(gB[-1]/gB[0]):.3f}",
    )

    # KEY
    print("\nKEY_NUMBERS (dual-run, rel tol 5%):")
    keys = {
        "maxRe_coupled": max_re,
        "slope_Omega": slope_A,
        "slope_Theta_to_G": slope_B,
        "slope_coupled": slope_C,
        "gainO_z0m2": float(gA[0]),
        "gainO_z0m10": float(gA[-1]),
        "gainThG_z0m2": float(gB[0]),
        "gainThG_z0m10": float(gB[-1]),
        "ratio_ThG_over_O_deep": ratio_deep,
    }
    for k, v in keys.items():
        print(f"  {k} = {v:.8g}")
    print(
        f"\nKEY maxRe={max_re:.6g} slopeO={slope_A:.6g} slopeThG={slope_B:.6g} "
        f"gO_m2={gA[0]:.6g} gO_m10={gA[-1]:.6g} "
        f"gThG_m2={gB[0]:.6g} gThG_m10={gB[-1]:.6g} ratio_deep={ratio_deep:.6g}"
    )

    print()
    if FAILURES:
        print(f"GATE FAIL: {FAILURES}")
        return 1
    print("GATE PASS: coupled resolvent scan (Ω / Θ→G / coupled gains)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
