---
id: OBS-01
title: With a bounded-at-edge anisotropic weight φ'(z)=2γ·½(1+tanh((z−z₀)/s)) (γ=0.4,z₀=−6,s=2), the BULK energy abscissa λ_max(−diag(φ')+sym_W(½Q_op)) is negative once a thin edge layer z<z_cut is excluded: −0.019/−0.106/−0.41/−0.65 at z_cut=−10/−8/−6/−4.
category: OBS
status: measured
confidence: 0.3   # lowered 2026-11-22: truncated-operator artefact, see CORR in body
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

**CORR 2026-11-22.** The numbers in this entry were produced with the discrete edge average truncated at `z_min`
(`build_blocks_z(edge_tail=False)`). With the exact edge tail (default since 2026-11-22) this entry's script prints
very different values (the large positive abscissas collapse), and a second review found the remaining values are
grid-specific because the discrete Hilbert transform is inaccurate near `z_min`. Treat the quantitative content of this
entry as a property of the truncated matrix, not of the operator. Status lowered from `verified` accordingly.
Details: `docs/findings_session_2026_11_22_qop_edge_tail.md`. The entry text below is left as written.

The numerical de-risking of (★) (LEM-06): the nonlocal forcing IS dominated by the
anisotropic transport dissipation in the bulk; what remains is the matched edge-layer
estimate. Evidence for `conjectured` LEM-06, does not itself close it.
