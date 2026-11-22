---
id: LEM-10
title: On L²_σ = L²(e^{−2σz}dz), σ∈(0,½), the constrained linearization M_σ = −2∂_z + K (K = e^{−z/2}T_c, c(s) = −tan(πs)/(2s+1), maximal realization) satisfies M_σ = 𝒩(−2∂_z − e^{−z/2})𝒩⁻¹, where 𝒩 is the Mellin multiplier N(s) = cos(πs)Γ(2s+2)/π, an isomorphism H^p_σ → L²_σ with p = 2σ+3/2. Consequently M_σ generates the explicit C₀-semigroup S(τ) = 𝒩T(τ)𝒩⁻¹, (T(τ)k)(z) = k(z−2τ)exp(−(1−e^{−τ})e^{−(z−2τ)/2}), with ‖S(τ)‖ = C*(σ)e^{−2στ} for every τ>0 (1 ≤ C*(σ) < ∞), and ‖(λ−M_σ)⁻¹‖ ≤ C(σ)/(Re λ+2σ) on the whole half-plane Re λ > −2σ.
category: LEM
status: proved-on-paper
confidence: 0.8
dependencies: [LEM-07, LEM-09]
supports: [LEM-06, THR-01]
cross_citations:
  - {label: "symbol identity", ref: "N(s+½)c(s+½) = −N(s)", note: "Lean: ConjugationSymbol.nSym_shift_mul_cSym"}
  - {label: "domain condition", ref: "N(½) = 0 ⟹ ĝ(½) = 0 on D(M_σ)", note: "Lean: ConjugationSymbol.nSym_half; recovers LEM-09(D) and explains the window σ < ½"}
  - {label: "replaces", ref: "docs/desk/w3_gearhart_pruss.md (withdrawn 2026-11-21)", note: "no Laplace inversion, no Gearhart–Prüss"}
verification:
  method: proof
  rounds: 1
  adversarial_vectors:
    - "independent Claude reviewer (fresh context, told to break it), 2026-11-21: steps 1–5 HOLD by own derivation; the domain identification D(M_σ) = 𝒩D(A) holds because the lines Re s = σ and σ−½ lie in a closed strip with no zero or pole of N (poles of 1/N at ±½ are never crossed)"
    - "end-to-end against the repo's real-space discretization (p.v. Hilbert quadrature, no Mellin input): closed-form evolution vs expm(τM_GG) 4.8e-4 at N=1500, 2.4e-4 at N=3000; transport alone 6.7e-2 (WC-3)"
    - "resolvent 𝒩(λ−A)⁻¹𝒩⁻¹ vs R_slit at six parameter points incl. −2σ < Re λ < 0 and complex λ: 1e-15…1e-29, except 1.6e-7 and 1.2e-9 at Re λ+2σ = 0.05, 0.1 (reviewer attributes to quadrature endpoint; not re-run regularised)"
    - "operator-norm probe (SVD over spans of u^a e^{−bu}, each member checked against the ladder equation): lower bounds C* = 1.0000 for σ ≤ 0.26, 1.0018 at 0.30, 1.1111 at 0.40 (converged); 1.39 at 0.45 and 2.41 at 0.499 (not converged, rising). No upper bound computed"
    - "reviewer CORRECTED the implementer: the OBS-09 modes map to k_λ ∈ H^p_σ with (λ−A)k_λ = πΓ(1+λ)𝟙 (Lemma D's defect; constant verified 2.81948734633 at λ=0.3), not to the homogeneous solution"
    - "strong continuity off the domain is only Hölder: ‖S(τ)g−g‖ ~ τ^{1−2σ} for ĝ(½) ≠ 0 (measured exponents 0.399–0.400 at σ=0.3, 0.099–0.100 at σ=0.45)"
    - "second independent Claude referee on the WRITTEN proof (manuscript Appendix app:conj), 2026-11-21: no step wrong; four small gaps in the domain lemma and five wording items, all repaired; defect constant re-verified to 12 digits"
    - "second referee did NOT reproduce the C* lower bounds above 1.06 (smaller family); the C* table is single-source"
    - "energy form: Re c(x+it) < 0 for all t when x in (0,1/2) (closed formula, re-derived by hand), so for sigma < 1/4 the form is non-positive on strip-regular data — indicates C* = 1 there; core argument not written"
    - "third independent Claude referee on Proposition prop:energy_form (written proof), 2026-11-21: could not refute; no finding above S3. Gaps repaired: the upper bound on N was used outside the lemma's stated range (now stated for compacts of (-1,inf)); formula for Re c undefined at x+it in 1/2+Z; two uncited semigroup facts; 'takes both signs' now proved. Its own checks: identity on two lines <= 5e-14 in 64 cases incl. sigma=0.03, 0.25; real-space finite-difference of e^{4 sigma tau}||S(tau)G||^2 against the form (sign and size agree, e.g. +1.6372e-2 vs +1.6381e-2 at sigma=0.4); generalized-eigenvalue search over five families at sigma in {0.02,0.1,0.2,0.25} found max ratio <= 1.000000000 and max form < 0"
    - "not confirmed by that referee: Engel-Nagel item numbers (II.1.3, II.1.7, II.1.10) against the book; it did not run the full probe"
  independent_rederivation: true
  machine_verification: "validate/edge_w3_conjugation_check.py (WC-1..WC-3, exit 0, ~30 s); reviewer scripts validate/review_2026_11_21_w3/ (not re-run by the implementer)"
lean_ref: "lean/SingularBlowup/ConjugationSymbol.lean — symbol identity, N(½)=0, semigroup law and characteristic ODE for the explicit flow, pointwise damping bound, sign of c on (0,½) and (½,1). NO operator statement is formalized."
tags: [track2, convective-weight, semigroup, W3, LEM-06-attack, review-2026-11-21]
---

Full argument and the reviewer's numbers: `docs/findings_session_2026_11_21_w3_conjugation.md`.

Why `proved-on-paper` and not `verified`: the Lean file checks only identities about explicit
functions. The analytic steps are on paper: the two-sided symbol bound
`|N(σ+it)| ≍ (1+|t|)^{2σ+3/2}`, the isomorphism `𝒩 : H^p_σ → L²_σ`, `D(M_σ) = 𝒩D(A)`, the `H^p`
multiplier bound for `Φ(z) = exp(−e^{−z/2})` (interpolation between `H⁰` and `H³`), strong
continuity, and the identification of the generator.

Consequences recorded here, each at the tier stated:

- **Spectrum** (derivation only; no numerical test left of the line): the closed half-plane
  `Re λ ≤ −2σ`, residual for `Re λ < −2σ` (`E = e^{λz/2}exp(−e^{−z/2})` annihilates the range of
  `λ − A`), no eigenvalues. The sentence "σ(M|L²_σ) = the line" in older notes is wrong for the
  maximal realization.
- **Not quasi-contractive for σ > ¼** (numerical lower bounds): `‖S(τ)‖` jumps to `C* > 1` at
  `τ = 0⁺`. Matches the instant weighted transient of review item E1.
- **Energy form (manuscript Proposition `prop:energy_form`; proved on paper).** On `D(M_σ)`:
  `Re⟨M_σG,G⟩_σ ≤ −2σ‖G‖²` iff `σ ≤ ¼`. So `C*(σ) = 1` and `‖S(τ)‖ = e^{−2στ}` for `σ ≤ ¼`, and
  `C*(σ) > 1` (not quasi-contractive) for `σ > ¼`. Probe WC-4, WC-5; Lean `reC_numerator_pos`,
  `reC_denominator_pos` for the sign factors only.
- **LEM-09 (C) needs a range restriction.** On the real axis `c < 0` on `(0,½)` and `c > 0` on
  `(½,1)` (Lean: `cRe_neg`, `cRe_pos`), so `Re c(σ+¼) > 0` at `t = 0` only for `σ > ¼`; LEM-09's
  sampled values all have `σ ≥ 0.27`. The reviewer reports `sup_t Re c(σ+¼+it) < 0` for
  `σ = 0.01…0.24` (grid `|t| ≤ 200`; not re-checked). The energy no-go therefore stands for
  `σ > ¼` only; whether `M_σ + 2σ` is dissipative for `σ < ¼` was not checked on the full domain.

Does not cover: the Gauge-(a) operator `−2∂_z + ½Q̃_op` (OBS-09), whose symbol acts on
`Re s ∈ (½, 3/2)` — conjugating it needs `1/N` across its pole at `s = ½`; and review item N1 (no
nontrivial solution is shown to satisfy the original system across `X = 1`).

Open: an upper bound on `C*(σ)` (reviewer's heuristic `(½−σ)^{−1/2}`, not fitted); equality with
`R_slit` on all of `L²_σ` as a proof; whether the manuscript's `M_σ` is literally the maximal
realization.
