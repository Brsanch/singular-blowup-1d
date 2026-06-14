---
id: LEM-01
title: The steady residual of eq. (2.4) vanishes identically on Y>1 for the singular profile — R_Ω = ā·Ω̄' − c_ω·Ω̄ − Θ̄' = 0 and R_Θ = ā·Θ̄' − (c_l+2c_ω)·Θ̄ = 0, with ā=2(Y-1), Ω̄=(Y-1)^{-1/2}, Θ̄=π/2.
category: LEM
status: verified
confidence: 1.0
dependencies: []
supports: [THR-01]
cross_citations:
  - {label: "Chen–Huang–Li 2026", ref: "arXiv:2604.01868 Thm 2.3", note: "the singular steady state we certify"}
verification:
  method: Lean
  rounds: 1
  adversarial_vectors: ["axiom audit (#print axioms): only propext, Classical.choice, Quot.sound — no sorryAx"]
  independent_rederivation: true
  machine_verification: "Lean: SingularBlowup.steady_residual_Ω / steady_residual_Θ; validate/verify_singular_steady.py (~2e-16)"
lean_ref: "lean/SingularBlowup/Profile.lean : steady_residual_Ω, steady_residual_Θ, hasDerivAt_Ωbar"
tags: [track1, backbone]
---

The exact (rather than ~2e-16 numerical, OBS-02) form of the steady-state gate.
`hasDerivAt_Ωbar` first establishes Ω̄'(Y) = −(1/2)(Y-1)^{-3/2} (rpow derivative).
The Ω-residual collapse uses (Y-1)·(Y-1)^{-3/2} = (Y-1)^{-1/2}; the Θ-residual uses
c_l+2c_ω = 0 (LEM-02) and Θ̄'=0. Zero sorry, zero axiom beyond Lean's standard three.
