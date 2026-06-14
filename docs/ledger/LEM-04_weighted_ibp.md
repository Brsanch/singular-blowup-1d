---
id: LEM-04
title: Weighted transport-dissipation identity — ∫_a^b (e^{-φ}G(-2G') + φ'e^{-φ}G²) dz = e^{-φ(a)}G(a)² − e^{-φ(b)}G(b)², i.e. ⟨G,-2∂_z G⟩_W = −∫φ'WG² − [WG²]_a^b.
category: LEM
status: verified
confidence: 1.0
dependencies: []
supports: [THR-01]
cross_citations:
  - {label: "Chen 2026", ref: "arXiv:2605.15130 §3", note: "anisotropic weighted estimates template; the φ' local rate"}
verification:
  method: Lean
  rounds: 1
  adversarial_vectors: ["axiom audit clean", "integrability via explicit ContinuousOn hypothesis on the integrand (honest side-condition, not assumed away)"]
  independent_rederivation: true
  machine_verification: "Lean: SingularBlowup.weighted_ibp"
lean_ref: "lean/SingularBlowup/Energy.lean : weighted_ibp"
tags: [track1, backbone]
---

The exact integration by parts behind the transport-dissipation step
(rigor_energy_estimate_setup.md §3): the local dissipation rate is φ'(z) and the
boundary term carries the convective gap (−2γ for φ=2γz), with the inflow piece
vanishing under the edge/Robin condition G(−∞)=0. Proved via FTC on F=−e^{-φ}G²
(intervalIntegral.integral_eq_sub_of_hasDerivAt). General φ, G, φ', G' (pointwise
derivatives + integrand continuity). This is the "needs no further work" term of
the energy estimate; the leap is the *forcing* term (LEM-06), not this. Zero sorry,
zero axiom beyond the standard three.
