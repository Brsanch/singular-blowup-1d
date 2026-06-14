---
id: OBS-04
title: Spectral-side continuation of M_GG — max Re eig(M_GG) ≈ −1e-3 with 0 eigenvalues at Re>0 across z_min∈[−18,−12], N∈[1400,2800], and → 0 as N→∞; this maximum is the MARGINAL transport essential spectrum on Re=0 (κ-stable under an edge g-mass damping that crushes the +40 field-of-values abscissa; delocalized, conjugate-paired), NOT a spectral gap and NOT the inadmissible (Y-1)^{-1/2} edge branch.
category: OBS
status: verified
confidence: 0.9
dependencies: []
supports: [LEM-06, THR-01]
cross_citations:
  - {label: "cap_scaffold.md", ref: "this repo", note: "essential spectrum sits on Re=0 unweighted; the gap is only after the weight"}
verification:
  method: numerics
  rounds: 1
  adversarial_vectors:
    - "weight-invariance / hi-k smoothness selection of 'the rightmost point eigenvalue' is FRAGILE on this non-normal operator (flips −0.7↔−22 with weighted-vs-unweighted scoring) — so NO specific point-eigenvalue value is claimed"
    - "edge-decay rate and c1=H(g)(1) readouts are unreliable (oscillatory modes; discrete-Hilbert edge bias OBS-03) — the membership B=0 question is NOT settled numerically"
    - "the near-0 cluster is κ-STABLE (does not move under −κ·diag(e^{−z}) even at κ=1e-2 that crushes the +40 abscissa) ⟹ delocalized essential spectrum, distinct from the inadmissible edge branch"
  independent_rederivation: true
  machine_verification: "validate/edge_spectrum_continuation.py"
lean_ref: none
tags: [track2-evidence, no-self-certification]
---

The discriminating run requested as the cheap next step. Two robust facts, one honest null:
(1) **no eigenvalue with Re>0** appears anywhere in the swept window — positive evidence for
stability; (2) but that maximum (`≈−1e-3`, the repo-wide "benign spectrum") is the **marginal
transport essential spectrum on Re=0**, NOT a gap — so the asymptotic-stability gap is
irreducibly a weighted-space property, and a direct unweighted eigenvalue computation cannot
exhibit it (it shows essential spectrum at Re=0). (NULL) the membership/Fredholm question
`B=H(g)(1)/(2λ)=0` is NOT settled numerically — the eigenvalue method inherits the same
edge-artifact / discrete-Hilbert ambiguity (OBS-03), which is itself confirmation that the
analytic Mellin/resolvent leap (LEM-06, `desk/edge_layer_estimate.md`) cannot be substituted by
numerics. **Repo-wide clarification:** the quoted "true spectrum ≈ −1e-3, stable" is the
*marginal essential spectrum*, not a spectral gap.
