---
id: LEM-07
title: The half-line (Y>1) Hilbert transform, in z=log(Y-1), is a convolution with Mellin symbol 𝓗(p)=cot(πp) (zero at p=−1/2 = the profile identity H[(Y-1)^{-1/2}]=0); the forcing ½Q_op is a −1/2 Mellin ladder with symbol c(ζ)=−tan(πζ)/(2ζ+1); hence M_GG G=λG is the first-order Mellin difference equation (2s+λ)a(s)=c(s+1/2)a(s+1/2).
category: LEM
status: proved-on-paper
confidence: 0.9
dependencies: [LEM-05]
supports: [LEM-06, THR-01]
cross_citations:
  - {label: "Gradshteyn–Ryzhik", ref: "3.222.2", note: "pv∫_0^∞ x^{a-1}/(1-x)dx = π cot(πa); H[(Y-1)^{-1/2}]=0"}
  - {label: "Erdélyi/Tricomi", ref: "one-sided (airfoil) Hilbert of powers", note: "(1/π)pv∫_0^∞ t^p/(x-t)dt = cot(πp)x^p — literature cross-check"}
verification:
  method: proof
  rounds: 1
  adversarial_vectors:
    - "4 independent re-derivations (keyhole contour; folding s→1/s quadrature, ~1e-12; Fourier-axis m(iξ)=−i coth(πξ), Re≡0; Erdélyi/Tricomi literature) — unanimous agrees_cot=yes (8/8)"
    - "strip −1<p<0 shown to be a TRUE convergence boundary (truncated integral diverges outside; out-of-strip value is analytic continuation, NOT a numerical artifact)"
    - "ζ=−1/2 (λ=+1) is an ORDER-2 pole (both tan and 1/(2ζ+1) blow up) = the time-translation symmetry mode in 𝒩 — sharpening, not a defect"
    - "forcing coefficient reproduced EXACTLY by sympy composition tracing build_blocks_z (difference simplifies to 0)"
  independent_rederivation: true
  machine_verification: "validate/edge_mellin_symbol.py (discrete-Hilbert numerics, <2% in-strip) + mpmath/sympy in workflow wf_0db34a7c (analytic, ~1e-10); NOT Lean"
lean_ref: none
tags: [track2, mellin, wiener-hopf]
---

The concrete analytic object the `edge_layer_estimate.md` rescope called for. Genuine
classical content (an explicit symbol computed from the operator), independently verified
by 4 methods — NOT a paraphrase. It reduces the strongly-non-normal nonlocal eigenvalue
problem to an explicit, single-valued, meromorphic **Mellin difference equation**. Full
derivation + the regime map + the circularity firewall: `docs/desk/edge_mellin_symbol.md`.
The symbol is `proved-on-paper` (4-method + mpmath cross-checked); the downstream *closure*
of stability (the Wiener–Hopf index of the difference equation) is the open LEM-06, NOT this.
