---
id: LEM-02
title: Edge / Θ-balance — Ū(1)=−2 (continuous across the edge), c_l+2c_ω=0, and the steady drift ā(1)=Ū(1)+c_l=0, so the singular term ā·Θ̄' has zero jump coefficient.
category: LEM
status: verified
confidence: 1.0
dependencies: []
supports: [THR-01, LEM-01]
cross_citations:
  - {label: "Chen–Huang–Li 2026", ref: "arXiv:2604.01868 eq. (3) / normalization", note: "c_l=−U(1), c_ω=−c_l/2 dynamic-rescaling pins"}
verification:
  method: Lean
  rounds: 1
  adversarial_vectors: ["axiom audit clean", "left-branch limit 2√(1-Y)-2 → -2 checked separately (continuity across the jump)"]
  independent_rederivation: true
  machine_verification: "Lean: SingularBlowup.Ubar_one / theta_balance / drift_edge_zero / theta_jump_coeff_zero"
lean_ref: "lean/SingularBlowup/Profile.lean : Ubar_one, Ubar_left_limit_at_one, theta_balance, drift_edge_zero, theta_jump_coeff_zero"
tags: [track1, backbone]
---

The consistency that lets the singular profile solve the *distributional*
Θ-equation across the X=1 jump: Θ̄' = (π/2)δ(X-1), and its coefficient
(Ū(1)+c_l)·(π/2) = 0 because Ū(1) = −2 = −c_l. Both branches of Ū agree at the
edge (√0=0), so Ū is continuous. Zero sorry, zero axiom beyond the standard three.
