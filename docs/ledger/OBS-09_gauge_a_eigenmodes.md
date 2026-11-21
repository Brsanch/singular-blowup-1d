---
id: OBS-09
title: For every Re λ > −2σ, σ ∈ (0,½), the Gauge-(a) operator −2∂_z + ½Q̃_op WITHOUT the domain condition Ĝ(½)=0 has the explicit eigenfunction G_λ = 𝓜⁻¹[Γ(2s+λ)Γ(2s+2)Γ(s)Γ(1−s)] ∈ L²_σ(ℝ) (contour max(0,−Re λ/2) < s0 < 1); G_λ is C^∞ up to the edge with G = O(v), is not edge-analytic (r_m = (−1)^{m+1}Γ(2m+λ)Γ(2m+2)), and has H[g](1) = −2Γ(1+λ) ≠ 0.
category: OBS
status: measured
confidence: 0.9
dependencies: [LEM-07]
supports: []
verification:
  method: numerics
  rounds: 2
  adversarial_vectors:
    - "real-space check with NO Mellin identity on the right-hand side (principal-value quadrature of the Hilbert transform): eigen-identity residual ≤ 2.2e-8 at λ = 0.3, 1, 3 (EM-1) and 1.6e-8 at λ = −1/4 (EM-5)"
    - "Codex review 2026-11-21: the first version allowed any contour 0<s0<1 and claimed G_λ → Γ(λ) for all Re λ > −2σ; both wrong for Re λ ≤ 0 (pole of Γ(2s+λ) at s = −λ/2). Reproduced at λ = −1/4 (A = 4.6506249067, G(16) = 29.4658775645) and corrected: far field is A_λ e^{−λz/2} + Γ(λ) + o(1) there"
    - "edge germ: r_1, r_2 match the residues at s = 1, 2 to 0.2%; |r_{m+1}/r_m| ~ (2m)^4, so zero radius of convergence — no contradiction with the edge-analytic nonlinear theorem (EM-3)"
    - "the repo's truncated-grid Gauge-(a) matrix does NOT carry the modes: no unstable eigenvalue, and H[g](1) of G_1 decays like e^{−T} where the continuum mode grows like e^{+T} (EM-4)"
    - "not tested in the probe: complex λ (Codex reports an independent complex-parameter residue check on commit 7514573)"
  independent_rederivation: false
  machine_verification: "validate/edge_gauge_a_eigenmode_check.py (EM-1..EM-5, exit 0, ~15 s)"
lean_ref: none
tags: [track2, gauge-a, eigenmodes, scope-of-THR-01, review-2026-11-21]
---

Context: found in an independent review of the manuscript on 2026-11-21; full statement, derivation and
the list of open questions are in `docs/findings_session_2026_11_21_gauge_a_eigenmodes.md`.

Why it follows: the ladder multiplier satisfies `(λ+2s)X(s;λ) = c(s+½)X(s+½;λ)` (LEM-07 ladder), `X` is
holomorphic and exponentially decaying on `max(0,−Re λ/2) < Re s < 1`, and
`𝓜[½Q̃_op G](s) = c(s+½)Ĝ(s+½)` holds without assuming `H[g](1) = 0`.

What it changes:

- **Point-spectrum uniqueness (LEM-08, W2) is a statement about the constrained domain.** The uniqueness
  argument kills the constant `P ≡ c_0` through `P(0) = 0`, which uses
  `(2∂_z+λ)G = e^{−z/2}T_cG = ½Q_op G`. For `½Q̃_op G = ½Q_op G − ½H[g](1)` the constant term has a Laplace
  pole at `s = 0`, so `P(0) ≠ 0` is allowed and `Ĝ = c_0 X` is a genuine eigenfunction.
- **Relation to the README firewall.** "`H(g)(1) = 0` is automatic for `L²_W` modes" is true for the
  frozen-exponent operator `−2∂_z + ½Q_op` (then `Q_op G → H[g](1)` at the edge, which is in `L²_σ`, `σ > 0`,
  only if it vanishes). It is NOT automatic for the Gauge-(a) operator with `Q̃_op`, which is the one the
  nonlinear theorem evolves.
- **Scope of the nonlinear theorem.** Edge analyticity excludes `G_λ`. In any finite-order `C^k` edge class
  the Gauge-(a) linearization has the growing solutions `e^{λτ}G_λ`, `Re λ > 0`.
- **Numerical evidence.** `#{Re λ > 0} = 0` on a grid truncated at `z_min` is evidence about the truncated
  problem: grid functions vanish near the edge, so the edge term in `d/dτ H[g](1)` (`−c_1/3`) is absent.

Not settled: whether the modes are physically relevant or an artifact of fixing the edge amplitude through
the unbounded functional `H[g](1)` (a different gauge might remove them).
