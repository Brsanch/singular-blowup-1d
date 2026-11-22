---
id: OBS-03
title: The +40 edge energy-abscissa is a continuum non-normal resonance carried by the g-unbounded (Y-1)^{-1/2} indicial branch, not a grid artifact and not a coercivity obstruction on the admissible class — A_full·e^{z_min/2}=0.099 const (z_min=−10..−16); bulk abscissa −0.41/−0.65 depth-independent; an infinitesimal g-mass penalty κ=1e-3 crushes +40 to negative; discrete Hilbert mis-resolves the profile at the edge node (−0.445·e^{−z_min/2} vs continuum 0).
category: OBS
status: measured
confidence: 0.3   # lowered 2026-11-22: truncated-operator artefact, see CORR in body
dependencies: []
supports: [LEM-06]
cross_citations:
  - {label: "Gradshteyn–Ryzhik", ref: "G-R 3.222.2", note: "H((Y-1)^{-1/2})=0 on the half-line — the continuum identity the discrete operator violates"}
verification:
  method: numerics
  rounds: 2
  adversarial_vectors:
    - "self-similar collapse A_full·e^{z_min/2}=0.099 is a CONTINUUM (Y-1)^{-1/2} rate, not a grid rate (would scale with N) — confirmed across 4 depths"
    - "N-drift at fixed z_min only ~15%/3× (sub-log) — the discrete-Hilbert echo riding on the continuum 0.099, not the leading term"
    - "g-penalty κ=1e-3 → A_full=−0.24 (101% crush): isolates the resonance to the g-unbounded branch"
  independent_rederivation: true
  machine_verification: "validate/edge_layer_obstruction.py (independent reproduction of the 6-lens workflow's decisive claims)"
lean_ref: none
tags: [track2-evidence, artifact-vs-real-settled]
---

**CORR 2026-11-22.** The numbers in this entry were produced with the discrete edge average truncated at `z_min`
(`build_blocks_z(edge_tail=False)`). With the exact edge tail (default since 2026-11-22) this entry's script prints
very different values (the large positive abscissas collapse), and a second review found the remaining values are
grid-specific because the discrete Hilbert transform is inaccurate near `z_min`. Treat the quantitative content of this
entry as a property of the truncated matrix, not of the operator. Status lowered from `verified` accordingly.
Details: `docs/findings_session_2026_11_22_qop_edge_tail.md`. The entry text below is left as written.

Settles the artifact-vs-real dispute that split the derivation lenses: BOTH a genuine
continuum profile-resonance AND a discrete-Hilbert contamination exist at the edge (same
e^{−z/2} rate, "two faces of the same unresolved profile singularity"). The decisive
disambiguators — depth-stable negative bulk and the κ-penalty crush — show the divergence
is the inadmissible (Y-1)^{-1/2} branch, excluded by admissibility. On the admissible
class the abscissa is negative but MARGINAL (≈−1e-3 = true spectral max), not coercive —
so the energy route gives no usable δ at the edge even with admissibility imposed. This is
the evidence base for the LEM-06 rescope (energy route dead; leap = non-energy spectral
bound). See `docs/desk/edge_layer_estimate.md`.
