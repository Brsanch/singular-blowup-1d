---
id: LEM-06
title: (★) The nonlocal forcing is dominated by the anisotropic transport dissipation — ⟨G, ½Q_op G⟩_W ≤ ∫φ'WG² − δ‖G‖²_W for G ⊥_W 𝒩 with G(−∞)=0, for some anisotropic weight W=e^{-φ} (φ'→2γ in the bulk) and δ>0.
category: LEM
status: conjectured
confidence: 0.60
dependencies: [LEM-04, OBS-01, OBS-03]
supports: [THR-01]
cross_citations:
  - {label: "Chen 2026", ref: "arXiv:2605.15130", note: "finite-codimension + anisotropic + IBP-twice template"}
verification:
  method: numerics
  rounds: 3
  adversarial_vectors:
    - "single exponential weight φ=2γz FAILS (field-of-values diverges, C_force 10→83) — refuted; the anisotropic weight is required"
    - "BULK abscissa < 0 and depth-stable (−0.41/−0.65 across z_min=−10..−16) — the bulk genuinely closes"
    - "6-lens derive + 18 adversarial verifies (2026-06-13): the EDGE layer does NOT close by ANY weighted-energy/Schur/flux method — the forcing inherits the (Y-1)^{-1/2} profile amplitude at the edge (A_full·e^{z_min/2}=0.099 const, reproduced); the +40 lives on the g-unbounded branch (κ=1e-3 penalty crushes it to negative)"
    - "circularity firewall: H(g)(1)=0 is AUTOMATIC for L²_W modes (not a new lemma); asserting it from the benign spectrum = the NS-C6/SQG circular move — declined"
  independent_rederivation: true
  machine_verification: "validate/energy_estimate_search.py (bulk<0) + validate/edge_layer_obstruction.py (self-similar collapse 0.099·e^{−z_min/2}; depth-stable bulk; g-penalty crush; discrete-Hilbert −0.445·e^{−z_min/2}); NOT a closed analytic proof"
lean_ref: none
tags: [track2, leap, critical-submission-blocker]
---

**The non-interpolative leap — RESCOPED 2026-06-13 (`docs/desk/edge_layer_estimate.md`).**
(★) is numerically established in the BULK (depth-stable coercive δ=0.02..0.78). The
edge-layer leg was attacked by 6 independent lenses + 18 adversarial verifications and
the verdict is a sharp NEGATIVE result: **no bounded-at-edge weighted-energy / Schur /
flux method closes the edge layer** — the nonlocal forcing inherits the base-profile
`(Y-1)^{-1/2}` amplitude there, and the only weight that would dominate it (φ'≳1 at the
edge) re-triggers the convective field-of-values divergence (cap_scaffold). The `+40`
edge abscissa is a continuum non-normal profile-resonance carried by the g-unbounded
`(Y-1)^{-1/2}` indicial branch (`Re λ∈[−1,0)`: L²_W-admissible but `g=e^{−z/2}G∉L²(dz)`),
excluded by admissibility — but only enforced numerically, not proven.

So the leap is NO LONGER an energy estimate: it is the **non-self-adjoint spectral /
Fredholm** statement "the L²_W realization of `M_GG ⊖ 𝒩` has no `Re λ ≥ 0` spectrum",
to be proven via a verified **resolvent bound** (cap_scaffold route 2) or the **Mellin /
Wiener–Hopf symbol** of `½Q_op` at the edge indicial exponent `r = −(λ+1)/2`, using the
continuum identity `H((Y-1)^{-1/2})=0` (which the discrete operator violates → numerics
mislead). Full route table, circularity firewall, and the two concrete next computations:
`docs/desk/edge_layer_estimate.md`. Stays `conjectured` (no-self-cert rule).
