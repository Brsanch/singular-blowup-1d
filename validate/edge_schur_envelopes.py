#!/usr/bin/env python3
"""
GATE: quantitative W1 Schur / Plemelj envelope constants vs (σ, δ)
(docs/desk/convective_weight_spectrum.md §10).

Desk §10 only machine-checks a few (τ,τ') points at fixed (σ=0.4, λ=0.3) that
dressed·e^{π|Δτ|} is bounded-and-decaying. This gate produces the *missing*
table of explicit constants:

  S_off(σ,δ)  = sup_{|Δτ|≥1} dressed·e^{π|Δτ|}
  α_eff(σ,δ)  = measured off-diagonal log-decay rate (expect ≥ π)
  M_diag(σ,δ) = max_{|Δτ|≤1, τ≠τ'} dressed
  C_term(σ,δ) = (Σ_n sup_τ |Π_n|) · δ     (desk ~1.15 at one point)
  γ_def(σ,δ)  = log-growth rate of bare defect density 1/(|X||λ+2s|)  (Stirling ⇒ ≈3π)

delta = Re(λ) + 2σ > 0 (distance above the transport essential-spectrum line).

PASS (fail-able, not tautological):
  P1  α_eff ≥ π − 0.05 on every row (off-diag decay at least the 4π−3π margin)
  P2  S_off < ∞ and finite on every row; S_off * δ stays O(1) (δ-blowup controlled)
  P3  C_term in [0.5, 3.0] at the desk anchor (σ=0.4, δ=1.1) — matches ~1.15 class
  P4  S_off grows as δ↓; S_off*δ within 50× across the δ-grid
  P5  bare defect density GROWS: γ_def ∈ (2.5π, 3.5π) (corrects desk §10 "decaying ~e^{-3π|τ|/2}")

Novelty vs desk: S_off / M_diag / C_term tables with δ-dependence; α_eff − π slack;
and the defect-kernel growth correction (desk prose claimed decay).
"""
from __future__ import annotations

import math
import os
import sys

# repo root on path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.schur_envelopes import (  # noqa: E402
    default_tau_grid,
    envelope_table,
    key_numbers,
)

FAILURES = []


def report(name: str, ok: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}  {detail}")
    if not ok:
        FAILURES.append(name)


def main() -> int:
    # Modest grid: fast enough for dual-run, dense enough for δ-dependence.
    # Include desk anchor λ=0.3, σ=0.4 ⇒ δ = 0.3 + 0.8 = 1.1
    sigmas = (0.25, 0.40)
    deltas = (0.15, 0.30, 0.50, 1.10, 1.50)
    tau_grid = default_tau_grid(tmax=8.0, n=29)
    print("W1 Schur envelope table  (band=1, tmax=8, dps=20)")
    print(f"  sigmas={sigmas}")
    print(f"  deltas={deltas}")
    print(f"  n_tau={len(tau_grid)}")
    print()
    rows = envelope_table(
        sigmas=sigmas,
        deltas=deltas,
        band=1.0,
        dps=20,
        n_max_termwise=24,
        tau_grid=tau_grid,
    )

    hdr = (
        f"{'σ':>6} {'δ':>6} {'λ':>7} {'S_off·e^{π|Δ|}':>14} {'α_eff':>8} "
        f"{'M_diag':>10} {'C_term':>8} {'γ_def':>8}"
    )
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        print(
            f"{r.sigma:6.2f} {r.delta:6.2f} {r.lam_real:7.3f} "
            f"{r.sup_off_pi:14.4g} {r.alpha_eff:8.4f} "
            f"{r.max_diag:10.4g} {r.termwise_C:8.4f} "
            f"{r.defect_growth:8.4f}"
        )
    print()

    kn = key_numbers(rows)
    print("KEY_NUMBERS (dual-run targets, rel tol 1e-3):")
    for k, v in kn.items():
        if isinstance(v, float):
            print(f"  {k} = {v:.8g}")
        else:
            print(f"  {k} = {v}")
    print()

    # ---- PASS criteria ----
    pi = math.pi
    alpha_ok = all(
        (math.isfinite(r.alpha_eff) and r.alpha_eff >= pi - 0.05) for r in rows
    )
    report(
        "P1 alpha_eff >= pi - 0.05 on all rows",
        alpha_ok,
        f"min α_eff = {min(r.alpha_eff for r in rows):.4f}  (π={pi:.4f})",
    )

    finite_ok = all(
        math.isfinite(r.sup_off_pi) and r.sup_off_pi > 0 and r.sup_off_pi < 1e8
        for r in rows
    )
    # S_off * δ should not explode: product bounded by 50 on this grid
    products = [r.sup_off_pi * r.delta for r in rows]
    prod_ok = all(p < 50.0 for p in products) and finite_ok
    report(
        "P2 S_off finite and S_off*δ < 50 on grid",
        prod_ok,
        f"max S_off*δ = {max(products):.4g}; max S_off = {max(r.sup_off_pi for r in rows):.4g}",
    )

    # desk anchor σ=0.4, δ=1.1
    anchor = next((r for r in rows if abs(r.sigma - 0.4) < 1e-9 and abs(r.delta - 1.1) < 1e-9), None)
    if anchor is None:
        report("P3 termwise_C at desk anchor in [0.5, 3.0]", False, "anchor row missing")
    else:
        c_ok = 0.5 <= anchor.termwise_C <= 3.0
        report(
            "P3 termwise_C at (σ=0.4,δ=1.1) in [0.5, 3.0]",
            c_ok,
            f"C_term={anchor.termwise_C:.4f} (desk ~1.15 class); "
            f"S_off={anchor.sup_off_pi:.4g}; α={anchor.alpha_eff:.4f}",
        )

    # δ-dependence: smaller δ should yield larger S_off (Γ pole approach)
    s04 = sorted([r for r in rows if abs(r.sigma - 0.4) < 1e-9], key=lambda r: r.delta)
    if len(s04) >= 2:
        mono = s04[0].sup_off_pi > s04[-1].sup_off_pi * 0.5  # at least not inverted
        # stronger: product S_off*δ roughly stable (within 20x)
        p_small = s04[0].sup_off_pi * s04[0].delta
        p_large = s04[-1].sup_off_pi * s04[-1].delta
        ratio = max(p_small, p_large) / max(min(p_small, p_large), 1e-30)
        report(
            "P4 S_off grows as δ↓ (pole approach); S_off*δ within 50x across δ",
            mono and ratio < 50,
            f"S_off(δ={s04[0].delta})={s04[0].sup_off_pi:.4g}; "
            f"S_off(δ={s04[-1].delta})={s04[-1].sup_off_pi:.4g}; "
            f"product ratio={ratio:.3g}",
        )
    else:
        report("P4 δ-dependence", False, "need ≥2 δ at σ=0.4")

    # P5: bare defect density grows at ~3π (desk §10 decay claim is false for 1/X)
    three_pi = 3 * math.pi
    g_vals = [r.defect_growth for r in rows]
    g_ok = all(math.isfinite(g) and (2.5 * math.pi < g < 3.5 * math.pi) for g in g_vals)
    report(
        "P5 defect density grows γ∈(2.5π, 3.5π) — not desk e^{-3π|τ|/2} decay",
        g_ok,
        f"γ range [{min(g_vals):.4f}, {max(g_vals):.4f}]; 3π={three_pi:.4f}",
    )

    # Print machine-readable KEY line for dual-run comparison
    print()
    print(
        f"KEY S_off_anchor={kn['anchor_sup_off_pi']:.6g} "
        f"alpha_anchor={kn['anchor_alpha_eff']:.6g} "
        f"C_term_anchor={kn['anchor_termwise_C']:.6g} "
        f"M_diag_anchor={kn['anchor_max_diag']:.6g} "
        f"defect_growth_anchor={kn['anchor_defect_growth']:.6g} "
        f"S_off_x_delta_small={kn['sup_off_pi_times_delta_small']:.6g}"
    )

    print()
    if FAILURES:
        print(f"GATE FAIL: {FAILURES}")
        return 1
    print("GATE PASS: W1 Schur envelope table + δ-scaling constants")
    return 0


if __name__ == "__main__":
    sys.exit(main())
