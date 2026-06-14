---
id: LEM-05
title: Edge indicial roots — (Y-1)^r solves the homogeneous edge ODE −2(Y-1)g'−g=λg iff r=−(λ+1)/2 (reduction −2(Y-1)g_r'−g_r = −(2r+1)(Y-1)^r); the weighted-admissible branch r≥−1/2 ⇔ λ≤0.
category: LEM
status: verified
confidence: 1.0
dependencies: []
supports: [THR-01]
cross_citations: []
verification:
  method: Lean
  rounds: 1
  adversarial_vectors: ["axiom audit clean", "iff proved both directions (not one-sided)", "real λ formalized; complex λ noted as the same computation on Re (un-formalized, documented honestly)"]
  independent_rederivation: true
  machine_verification: "Lean: SingularBlowup.indicial_reduction / indicial_root_iff / indicial_admissible"
lean_ref: "lean/SingularBlowup/Energy.lean : indicial_reduction, indicial_root_iff, indicial_admissible, hasDerivAt_gpow"
tags: [track1, backbone]
---

The rigorous form of the Robin/edge condition (findings_edge_condition_derivation.md
§5): the indicial polynomial −(2r+1) fixes r=−(λ+1)/2, and G=(Y-1)^{1/2}(Y-1)^r is
weighted-admissible (decays into the edge) iff r≥−1/2 iff Re λ≤0 — so the inadmissible
homogeneous branch is dropped exactly when Re λ>0, which is the 𝒩-complement /
finite-codimension structure THR-01 needs. Lean covers real λ; the complex case is the
identical computation on the real part. Zero sorry, zero axiom beyond the standard three.
