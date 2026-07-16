---
id: OBS-07
title: Depth / weak-coupling scan of the exact defect-closed ladder resolvent in L²_σ (σ=0.4) — gain(λ,z₀) = ‖(λ−M)⁻¹F‖_σ/‖F‖_σ for F = e^{−cosh(z−z₀)}, z₀ = −2…−14, λ ∈ {0.3, 0.05, −0.2, −0.5, 0.3+1.5i, −0.7} (transport line at −0.8). VERDICT — gain ∝ e^{z₀/2} = ε, λ-uniform to 0.2% incl. dist-0.1 from the line: the response to edge forcing VANISHES linearly in the weak coupling (measured 0.30/0.082/0.0134/0.00503/0.00187 at z₀ = −2/−6/−10/−12/−14; log-slopes −0.33, −0.45, −0.490, −0.494 → ½). Strong numerical support for σ(M|_{L²_σ}) = {Re λ = −2σ} with ‖R(λ)‖ ~ C/dist attained on bulk-scale forcing (where the termwise bound gives ≈1.15/dist).
category: OBS
status: measured
confidence: 0.85
dependencies: [LEM-09]
supports: [LEM-06]
verification:
  method: numerics
  rounds: 2
  adversarial_vectors:
    - "continuum-exact: Gamma/Mellin/Bessel-K arithmetic on the line only — no discrete Hilbert (the repo's proven edge-bias trap cannot enter); cancellations to 10^67 at dps 90 (~23-digit margin)"
    - "cross-implementation anchor: gain(0.3, −6) = 0.08173 reproduced independently by a float64 numpy implementation that also PDE-verifies its G against the continuum operator ((λ−M)G−F residual 1.4e-3)"
    - "CAUGHT + FIXED a dps-STABLE analytic bug: the defect constant 𝒢₀(½) evaluated by ±δ symmetric limit is garbage at depth (δ²·e^{(3/2)n*ln n*} curvature swamps the value; passed a 90/160/260-digit consistency test; produced spurious 1e12–1e15 gains at z₀ ≤ −12; the line-eq self-check is blind to the fix constant). Only the exact closed form Π_n(½) = X(½)/[X(½+n/2)(λ+1+n)] (even n) is safe"
    - "mechanism cross-check: the e^{z₀/2} law matches Lemma C's e^{|a|/2}-large dissipative K-form on edge bumps; the singular weak-coupling discontinuity lim‖R_εF₀‖ = 0 ≠ ‖R_AF₀‖ is explained by high X-shaped rungs dominating the ladder (C₀(line) ≈ C₀(½) to 15+ digits at depth)"
  independent_rederivation: true
  machine_verification: "validate/edge_convective_resolvent_scan.py (--deep reproduces this grid; quick mode = regression subset)"
lean_ref: none
tags: [track2, convective-weight, numerics, LEM-06-probe]
---

Full table + reading: `docs/desk/convective_weight_spectrum.md` §7. Interpretation contract
(fixed before the run): bounded dressing ⟹ supports the LEM-06 essential-spectrum estimate;
`e^{θ|z₀|}` growth ⟹ refutes. Outcome: not merely bounded — *decaying at the exact
weak-coupling rate*. NOT a proof: the analytic object remains the §8 bent-contour Barnes bound;
this measurement de-risks it (the hard part — deep-edge data — is where the true solution is
smallest, and the numerically-exhibited two-stage cancellation `𝒢 = X·[C₀(s) − C₀(½)]` is the
structure the proof must capture).
