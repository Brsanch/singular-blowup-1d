# Independent hostile review — `cf2f6b2` (Track A Lemma A3, analytic energy closure)

**Date:** 2026-11-22 · **Review form:** independent — single Claude session (Opus 5.5), memnet off; analytic re-derivation + numerical probes. Not the implementer. No tracked repo files modified.
**Gate:** `reslock run -- python3 validate/edge_track_a_proof_lemmas.py` → `ALL PASS (TA-1..TA-5)` as committed.

## Verdict
No S1. Lemma A3 plausibly closes the *a priori* energy bound (1.3), conditional on two repairs: the missing commutator in (2.9) (numerically shown harmless, not yet written as a proof) and the zero-margin γ. Theorem 2 remains an a priori estimate only; existence is unproved and the surrounding document still claims it.

## Findings

### S2 — (2.9) is not literally true: missing commutator
⟨D⟩ is a Fourier multiplier along the line; u = e^{-x/2}e^{∓iρ/2} shifts the Mellin line, so ⟨D⟩^{-1/2}u ≠ u⟨D⟩^{-1/2}. The correct equation is v_s + ⟨D⟩^{-1/2}u⟨D⟩^{1/2}v = …, not v_s + uv. The extra term has the same weight as uv and is only O(1/ξ) lower at high frequency, so it is not an absorbable remainder. (2.10) and `TA-5` treat the commutator-free operator.

**Repair (numerically checked):** after conjugation by e^{-σx} everything is constant-coefficient, so Re⟨Lv,v⟩_{H^k_hier(μ)} has the exact symbol
Re[e^{-iρ/2} · conj(K̂(ξ − i/4)) · Σ_m c_m (iξ+σ−½)^m conj(iξ+σ)^m],  K̂(η) = (1+η²)^{-1/2}.
Min ratio to ‖v‖²_{H^k(σ)} over |ρ| ≤ π/2, all ξ: **0.598** (k=2, B=6, μ=0), **0.514** (k=3, B=8, μ=3/8). The commutator adds at most 7.24° of phase. Coercivity survives; the proof needs this symbol bound written out.
Side effect: ⟨D⟩ is line-dependent, so "v" on μ=0 and on μ=3/8 are two different functions; ℰ(s) sums them (harmless, S3).

### S2 (scope) — document still claims results it does not prove
Theorem 2's own parenthetical scope note is honest (a priori bound for any smooth solution; existence/uniqueness deferred). But the title ("Global Nonlinear … Compatible Manifold"), the Status line ("Closed … proof"), and Lemma A4 ("reconstructing the global classical solution … ∎") assert existence, which does not follow. Construction is the hard part: the regularisation must preserve the nonlocal m = 0 and m = ½ cancellations.

### S3 — γ has zero margin
With γ = 12C_eq/(a c₁ c_u′), the remaining dissipation γc_u′/(2C_eq)·‖⟨D⟩^{-1/2}W_reg‖² equals 6/(a c₁) exactly — all of it consumed by the G-equation. Nothing is left for the rank-1 Young step or for the ∫‖⟨D⟩^{-1/2}W‖² term in (1.3). Fix: γ = 24C_eq/(a c₁ c_u′).

### S3 — pole cancellations: signs correct, one inference unjustified
With N_g(m) = Γ(2m+2)cos(πm)/π (recovered from Φ̂_* = Γ(2m+2), φ̂_* = π/cos(πm)): Res_{m=½}(F̂_A/N_g) = −(π/2)a((QG̃)W); rank-1 term gives +(π/2)a((QG̃)W); they cancel. πċ₀ cancels exactly. h_reg(·,0)=0 since a(W₀) = −2c₀ at s=0.
But "forcing pole-free ⇒ h_reg holomorphic" does not follow as stated: multiplication by u shifts the Mellin variable by ½, so a pole of ĥ_reg at ½ is coupled to ĥ_reg at m=1. Harmless: the energy runs in physical space and each forcing only needs the lines it is paired on (−¼, ⅛ for the source; 0, ⅜ for Term A and rank-1).

### S3 — product/multiplier estimates asserted, not proved
‖⟨D⟩^{½}R̃‖_{−¼} + ‖⟨D⟩^{¾}R̃‖_{⅛} ≤ C‖G̃‖‖⟨D⟩^{½}G̃‖; QG̃ bounded multiplier on H^{½}, H^{¾}; |a((QG̃)W)| ≤ C‖G̃‖‖⟨D⟩^{-1/2}W‖. Frequency-by-frequency counts work for every term, including Q(P(G̃)) (Q recovers the derivative D(QG̃)·G̃ loses). Term A derivative counts check: W_reg ∈ H^{½} on μ=0, H^{¾} on μ=⅜; pairings need v ∈ H^{k+½}, supplied by strip dissipation. Correction: D[a∘P](G̃)·W needs ‖⟨D⟩^{½}G̃‖, not (2.1); gives ‖β‖_{L¹} = O(ε²/√a), better than the stated O(ε²/a).

### S3 — `TA-5` is weaker than it looks
- Convexity argument correct: 2001-point sweep on [0, π/2] has its minimum exactly at the endpoint (0.4840, 0.3987).
- But it checks a pointwise *sufficient* matrix for the commutator-free operator, not "the exact Hermitian matrix of (2.10)".
- Threshold cos(π/4)/2 is unmotivated: k=3, B=6 gives λ_min = 0.3528 > 0 yet fails the gate.
- Insensitive past ρ = π/2 until ρ ≈ 2.5 (k=3) / 3.0 (k=2); never tests such ρ anyway. No mutation control.

### doc
- (2.8): should read 2C₁a, not C₁a.
- ‖v‖_{H^k(σ)} ≳ ‖⟨D⟩^{-1/2}W_reg‖ is an equivalence, not "exact": ⟨D⟩_μ on line σ has symbol (15/16 + ξ² − iξ/2)^{-1/4} ≍ ⟨ξ⟩^{-1/2}.
- e^{±x/4} shift fine up to constant c_shift (triangular invertible constant-coefficient map).

## Checked and correct
(2.8) cosh-weight Plancherel normalisation and |ξ|tanh(2ρ|ξ|) ≥ c₁⟨ξ⟩ − C₁; N_g order 2σ+3/2 (2 on ¼, 11/4+¼ = 3 on ⅝); N_r⁻¹ orders −1 (on −¼) and −7/4 (on ⅛); ‖N_g[e^{-ut}φ_*]‖ ~ t^{-1/4} on ⅝ and uniformly bounded on μ < ½, Young L²∗L¹→L² for W_rk1; W_alg = O(‖G̃‖²); bootstrap leftovers O(ε⁴/a²), O(ε⁴/a³); absorption needs ε ≲ a c₁.

## Not checked
Rigorous weighted strip-boundary product estimates; uniformity of C_eq, c_shift; existence of a smooth solution in this class (vacuity); Lemmas A1, A2, A4; Lean file; `NM-1..3`; that N_g is exactly this Mellin multiplier (only consistency with Φ_*, φ_* checked).

Probe script (temporary scratchpad, not in repo): `probe_ta5.py`, an exact-symbol coercivity check with/without the commutator plus TA-5 rho/B mutations.
