---
id: LEM-03
title: c_ω-absorption — in the edge-resolved variable g(Y)=(Y-1)^{-1/2}G(log(Y-1)), the transport −ā g_Y (ā=2(Y-1)) plus the source c_ω g (c_ω=−1) collapses to −2∂_z G: −2(Y-1)g'(Y)−g(Y) = −2(Y-1)^{-1/2}G'(log(Y-1)).
category: LEM
status: verified
confidence: 1.0
dependencies: []
supports: [THR-01]
cross_citations: []
verification:
  method: Lean
  rounds: 1
  adversarial_vectors: ["axiom audit clean", "product+chain rule derivative of (Y-1)^{-1/2}G(log(Y-1)) built and reduced via rpow_add + mul_inv_cancel"]
  independent_rederivation: true
  machine_verification: "Lean: SingularBlowup.cω_absorption; numerical 1e-11 (src/stability_edge_z.py)"
lean_ref: "lean/SingularBlowup/Profile.lean : cω_absorption, hasDerivAt_gOf, gOf"
tags: [track1, backbone]
---

The "+G/−G cancellation": the +G produced by −ā g_Y (from differentiating the
(Y-1)^{-1/2} factor) exactly cancels the −G from c_ω g, leaving the clean transport
generator −2∂_z. This is what makes M_GG = −2∂_z + ½Q_op (rigor_energy_estimate_setup.md
§1). Proved for arbitrary C¹ G (hypothesis: ∀ z, HasDerivAt G (G' z) z), not just the
profile. Zero sorry, zero axiom beyond the standard three.
