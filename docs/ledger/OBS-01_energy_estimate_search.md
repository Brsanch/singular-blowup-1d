---
id: OBS-01
title: With a bounded-at-edge anisotropic weight φ'(z)=2γ·½(1+tanh((z−z₀)/s)) (γ=0.4,z₀=−6,s=2), the BULK energy abscissa λ_max(−diag(φ')+sym_W(½Q_op)) is negative once a thin edge layer z<z_cut is excluded: −0.019/−0.106/−0.41/−0.65 at z_cut=−10/−8/−6/−4.
category: OBS
status: verified
confidence: 1.0
dependencies: []
supports: [LEM-06]
cross_citations: []
verification:
  method: numerics
  rounds: 1
  adversarial_vectors: ["global abscissa ≈+40 is a grid-scale edge artifact (worst mode at edge node, grows 35.7→41.9 with N=700→2000) — excluded by the Robin layer"]
  independent_rederivation: false
  machine_verification: "validate/energy_estimate_search.py"
lean_ref: none
tags: [track2-evidence]
---

The numerical de-risking of (★) (LEM-06): the nonlocal forcing IS dominated by the
anisotropic transport dissipation in the bulk; what remains is the matched edge-layer
estimate. Evidence for `conjectured` LEM-06, does not itself close it.
