---
id: OBS-02
title: The singular steady state (Thm 2.3) is independently certified numerically — steady residual of eq. (2.4) ~2e-16 on X>1 with (c_l,c_ω)=(2,−1); U(1)=−2 exact; ℝ-Hilbert tool validated vs Lorentzian to 1.4e-6.
category: OBS
status: verified
confidence: 1.0
dependencies: []
supports: [LEM-01, LEM-02]
cross_citations:
  - {label: "Chen–Huang–Li 2026", ref: "arXiv:2604.01868 Thm 2.3"}
verification:
  method: numerics
  rounds: 1
  adversarial_vectors: ["~15% error in the numeric Hilbert of the singular profile is far-field 1/√X truncation, not a steady-state failure (analytic-U residual is 2e-16)"]
  independent_rederivation: true
  machine_verification: "validate/verify_singular_steady.py (PASS)"
lean_ref: none
tags: [profile-evidence]
---

The numerical gate LEM-01 makes exact. Independent confirmation of a ~2-month-old
single-group claim (one of the two robust, hard-to-scoop contributions per README).
