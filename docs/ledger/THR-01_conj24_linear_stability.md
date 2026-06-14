---
id: THR-01
title: The singular Hou–Luo self-similar profile (Ω̄=(Y-1)^{-1/2}, Θ̄=(π/2)1{Y>1}, c_l=2, c_ω=-1) is asymptotically stable modulo the neutral subspace 𝒩={ξ@1, μ, ν} (Conjecture 2.4).
category: THR
status: speculative
confidence: 0.70
dependencies: [LEM-06]
supports: []
cross_citations:
  - {label: "Chen–Huang–Li 2026", ref: "arXiv:2604.01868 Conj. 2.4", note: "numerical-only; no rigorous stability spectrum published"}
verification:
  method: numerics
  rounds: 3
  adversarial_vectors:
    - "leading apparent unstable mode is a free-boundary edge artifact (δ^{-3/2}, grid-scale) — confirmed by N-convergence, not a converged eigenvalue"
    - "transient (non-normal) growth — ruled out via edge-conditioned operator, μ_an<0 resolution-independent"
    - "coupling-driven instability — coupled (G,h) shape march also relaxes"
  independent_rederivation: false
  machine_verification: none
lean_ref: none
tags: [headline, critical-submission-blocker]
---

**Claim.** The linearization 𝓛 of eq. (2.4) about the singular profile generates a
strictly contractive semigroup on H¹_W ⊖ 𝒩, i.e. the profile is asymptotically
stable modulo the dilation/amplitude/time-translation symmetries.

**Why speculative, not verified.** Per the no-self-certification rule, the repo's
central claim cannot exceed `speculative` until external review. The *numerical*
evidence is strong (no Re>0 mode at two discretizations; convective gap 2γ; no
transient growth; finite-amplitude + fully self-consistent co-evolving relaxation;
paper-matched normalization). But the proof reduces to LEM-06 (the (★) bulk
coercivity of the nonlocal forcing) + the H¹_W coupling closure, which are
`conjectured`, not proved. The Track-1 algebraic backbone (LEM-01..05) is
Lean-verified but is the *skeleton*, not the estimate — formalizing the backbone
does not raise this entry's status.

**Reduction (rigor_energy_estimate_setup.md §7).** 𝓛 contractive ⇐ (★) holds for
G ⊥_W 𝒩 with the edge condition G(−∞)=0, plus its ∂_z-level analogue. The transport
dissipation (LEM-04), the edge balance making the profile a steady solution
(LEM-01..02), the c_ω-absorption that gives the clean −2∂_z generator (LEM-03), and
the indicial admissibility that drops the inadmissible homogeneous edge branch
(LEM-05) are the closed pieces this reduction stands on.
