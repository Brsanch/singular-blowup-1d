---
id: OBS-08
title: High-frequency corner probe — for F = e^{−cosh(z−z₀)}e^{iΩz} the dressing gain/transport-gain rises monotonically to 1 with |Ω| at every depth (z₀=−2: 0.457→1.0000, z₀=−6: 0.123→0.9999, z₀=−10: 0.020→0.9985 across Ω = 0/−5/−15/−40; ±Ω symmetric): the region the slit-plane crude Schur could not close (and OBS-07 never probed) is benign — high-frequency edge forcing feels the weakened coupling c(σ+iΩ) ~ 1/(2|Ω|) and responds as pure transport. Nowhere in the probed (z₀, Ω) plane does the dressing exceed ~1.
category: OBS
status: measured
confidence: 0.85
dependencies: [LEM-09, OBS-07]
supports: [LEM-06]
verification:
  method: numerics
  rounds: 1
  adversarial_vectors:
    - "probes exactly the quadrant where the pre-subtraction Schur analysis failed (input |τ'| ≫ |τ|, same sign) — falsification-first before writing the closure"
    - "±Ω symmetry verified (conjugate structure); the low-frequency suppression crossover (z₀=−10, Ω=−5: dressing 0.44) matches the §3 band mechanism"
    - "same continuum-exact machinery as OBS-07 (Bessel-K/Gamma line arithmetic; closed-form defect fix; no discrete Hilbert)"
  independent_rederivation: false
  machine_verification: "validate/edge_slitplane_representation.py (P3 pins dressing(z₀=−6, Ω=−15) = 0.9968)"
lean_ref: none
tags: [track2, convective-weight, numerics, LEM-06-probe]
---

Context: the slit-plane Plemelj reduction (desk §10) closes the dressed-Schur bound in all
four (τ,τ')-quadrants once the split-regularized kernel is used; this probe was run BEFORE
that resolution was found, as the falsification test of the then-open corner. Both now agree:
the corner is benign (analytically via the subtracted kernel's e^{−4π(τ−τ')} decay,
numerically via dressing → 1). Full table: `docs/desk/convective_weight_spectrum.md` §10(c).
