---
id: LEM-06
title: (★) The nonlocal forcing is dominated by the anisotropic transport dissipation — ⟨G, ½Q_op G⟩_W ≤ ∫φ'WG² − δ‖G‖²_W for G ⊥_W 𝒩 with G(−∞)=0, for some anisotropic weight W=e^{-φ} (φ'→2γ in the bulk) and δ>0.
category: LEM
status: conjectured
confidence: 0.60
dependencies: [LEM-04]
supports: [THR-01]
cross_citations:
  - {label: "Chen 2026", ref: "arXiv:2605.15130", note: "finite-codimension + anisotropic + IBP-twice template"}
verification:
  method: numerics
  rounds: 2
  adversarial_vectors:
    - "single exponential weight φ=2γz FAILS (field-of-values diverges, C_force 10→83) — refuted; the anisotropic weight is required"
    - "global energy abscissa ≈+40 but the worst mode is the grid-scale edge singularity (grows with N) — not a converged eigenvalue"
  independent_rederivation: false
  machine_verification: "validate/energy_estimate_search.py (BULK abscissa < 0: −0.019/−0.106/−0.41/−0.65 for z_cut=−10/−8/−6/−4); NOT a closed analytic proof"
lean_ref: none
tags: [track2, leap, critical-submission-blocker]
---

**The non-interpolative leap.** This is the genuine remaining analytic work, NOT
formalized in Track 1. (★) is numerically established *in the bulk* (energy abscissa
< 0 with a bounded-at-edge anisotropic weight) but reduces to: (a) a rigorous
matched **edge-layer estimate** at Y→1 where G→0 (Robin) and the forcing is regular
(~½H(g)(1)(Y-1)^{-1/2}), matched to the bulk coercivity; and (b) the H¹_W coupling
closure (differentiate, control M_Gh and [∂_z,F]). Tool: an explicit Hilbert/averaging
**kernel estimate** on Q_op G = δU/(Y-1), δU=∫_1^Y H(g) — not the pessimistic
field-of-values. Stays `conjectured` no matter how many gates pass (no-self-cert rule).
Attack doc belongs in `docs/desk/` (OPEN-GAP convention).
