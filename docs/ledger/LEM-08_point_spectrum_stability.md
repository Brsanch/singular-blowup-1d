---
id: LEM-08
title: The frozen Ω-shape generator M_GG has NO Re λ≥0 point eigenvalue modulo 𝒩 — an exact integer obstruction. With Ĝ_p=Γ(s+λ)Γ(s+2)Γ(s/2)Γ(1−s/2) and the only homogeneous freedom a 1-periodic entire P (Laurent in w=e^{2πis}): line-L² forces degree d=0 (k_max<3/4), but killing the s=0 edge-constant pole (residue 2Γ(λ)) forces d≥1; incompatible ⟹ no point eigenfunction except at λ∈ℤ (= 𝒩, where the recurrence terminates to a polynomial).
category: LEM
status: proved-on-paper
confidence: 0.85
dependencies: [LEM-07]
supports: [THR-01]
cross_citations:
  - {label: "Liouville/Laurent rigidity", ref: "zero-free 1-periodic entire ⟹ constant", note: "no theta-ratio escape: the degree-vs-growth law is exact"}
verification:
  method: proof
  rounds: 2
  adversarial_vectors:
    - "construction machine-exact (~1e-49): τ=π/sin(πs/2)=Γ(s/2)Γ(1−s/2); Ĝ_p 4-Gamma; line-decay −3π/2 (diff-quotient); s=0 residue 2Γ(λ); WH split clean; eigenfunction = Fox-H H^{1,3}_{3,1} (bitwise identical)"
    - "3 independent routes (residue/pole-collision, Pincherle/Liouville, direct numerical L² scan) + the in-strip arg-principle count ALL converge on no Re λ≥0 point eigenvalue mod 𝒩"
    - "corrected TWO bookkeeping errors in the heuristic (far/edge pole assignment reversed; sin(πs) anti-periodic INADMISSIBLE so per-zero growth is 2π not π) — conclusion survives, sharper (margin holds at d=1)"
    - "firewall CLEAN: the const∉L²_W link is a property of the weight (Mellin-L² ⟺ no s=0 pole), computed from Ĝ_p, not injected from the benign spectrum (validate/edge_synthesis_firewall_audit.py)"
  independent_rederivation: true
  machine_verification: "validate/edge_pointspectrum_rigid.py (int|Ĝ_p·P|²=22.6 for d=0, divergent ∀d≥1) + edge_pole_collision.py + edge_route3_periodic_correct.py + edge_tau_independent_verify.py"
lean_ref: none
tags: [track2, point-spectrum, fox-h]
---

The payoff of the WH-factorization arc: with the tan-Barnes elementary (LEM-07/§5.3), the
point-spectrum question becomes the EXACT integer obstruction (§5.4): no 1-periodic entire `P`
can make `Ĝ_p·P` simultaneously line-`L²` (forces `d=0`) and `s=0`-pole-free (forces `d≥1`). So
**no `Re λ≥0` point eigenvalue mod `𝒩`** — the frozen Hou–Luo edge generator is **point-spectrally
stable**.

**This is NOT a strict spectral gap (so NOT yet Conjecture 2.4).** Three honest scope limits: (1) the
**essential spectrum on `Re=0`** (`R(0)=0 ∀λ` ⟹ non-Fredholm line) is untouched and marginal — a
strict gap `Re σ≤−δ` still needs the edge-unbounded weight that re-triggers `cap_scaffold` (LEM-06,
the genuine remaining object); (2) the `const∉L²_W` edge-mass step deserves a clean Lean/interval
statement (the now-*computed* content of the former firewall); (3) frozen generator — the
nonlinear/co-moving correction is a separate leg. Status `proved-on-paper` (machine-gated structural
argument), not `verified` (no Lean), and point-spectrum only. Full detail: `docs/desk/edge_mellin_symbol.md`
§5.3–5.5.
