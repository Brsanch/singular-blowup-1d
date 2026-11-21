/-
# LinearStrictGapClosure.lean — arithmetic bookkeeping for the linear resolvent argument (W1, W2)

**Scope (corrected 2026-11-21).** Every theorem in this file is a statement about real numbers
and integers. Operator-norm bounds enter as *hypotheses*; no Gamma function, singular integral,
weighted space, resolvent or semigroup is formalized here. The theorems record how externally
supplied bounds combine. They do not certify the linear stability result of the manuscript.

What changed since this file was written (2026-11-20):
  • **W3 is withdrawn.** Semigroup generation and the Gearhart–Prüss decay are no longer claimed
    (manuscript Remark 1; `docs/desk/w3_gearhart_pruss.md`, Correction). The two `w3_*` theorems
    below remain true arithmetic facts; their intended application does not go through.
  • **W1 is a vertical-line bound.** The analytic bounds hold on each line `Re λ = -2σ + δ`,
    uniformly on compact `δ`-intervals, not on whole half-planes `Re λ ≥ -2σ + δ`.
  • **W2 holds on the constrained domain only** (`Ĝ(1/2) = 0`). Without it the Gauge-(a) operator
    has an explicit eigenfunction for every `Re λ > -2σ`
    (`docs/findings_session_2026_11_21_gauge_a_eigenmodes.md`).

Declaration names are unchanged (they are cited elsewhere); names such as `certified_…` and
`…_master…` describe the intended use, not what Lean checks.

Original description (2026-11-20), kept for reference — read it with the scope note above.
It chains the three analytical walls of the linear strict-gap argument for the singular
1D Hou–Luo self-similar blow-up profile (`Ω̄(X) = 1_{X>1}/√(X-1)`, `c_l = 2`, `c_ω = -1`)
into arithmetic statements with in-repo rational bounds:

  • **W1 (`SlitSchurMargin` Route F):**
    Combines the 4-quadrant Stirling/split-Cauchy exponential margin
    (`lem_stir_schur_four_quadrant_margin`), the unitary-phase diagonal Riesz
    commutator bound (`modulus_commutator_cancels_cauchy_pole`), and the rank-one
    defect ω-cancellation (`slit_defect_rank_one_margin`) into a uniform vertical-line
    resolvent bound `M_tot(σ, δ) = M_diag / δ + M_far / δ + M_defect / (1 + |ω|)` on
    the line `Re λ = -2σ + δ` (the three component bounds are hypotheses here).

  • **W2 (`SlitSchurMargin` Left-Strip Bridge + `PointSpectrum` LEM-08):**
    Combines the left-strip Paley–Wiener embedding exponents
    (`left_strip_edge_exponent_le`, `left_strip_far_exponent_le`,
    `shift_half_weight_exponent_eq`), the pole-free Gamma real-part floors on
    `[σ - 1/2, σ]` (`left_strip_gamma_real_floors`), and the discrete integer
    Laurent obstruction (`PointSpectrum.no_l2_and_edge_kill`) to prove that any
    admissible left-strip solution has trivial periodic multiplier (`kMax = 0 ∧ kMin = 0`)
    and cannot support an edge-pole-killing Laurent mode.

  • **W3 (`SlitSchurMargin` Reference-Transport Perturbation & Frequency Pincer):**
    Uses `transport_coupling_frequency_pincer` (`|ω| ≤ |2τ + ω| + |2τ|`) to prove that
    every frequency `τ ∈ ℝ` lies either in the high-transport regime (`|2τ + ω| ≥ |ω|/2`)
    or the high-coupling regime (`|2τ| ≥ |ω|/2`). WITHDRAWN 2026-11-21: this was meant to
    yield `O((1 + |ω|)^{-1})` decay of `δR(λ) = R(λ) - R₀(λ)` and hence Gearhart–Prüss
    decay. The operator bound it needs (`L²_σ → L²_{σ+1/2}`) is not established, and a
    `1/(1+|ω|)` operator norm is not integrable in `ω`. Only the arithmetic is checked.

  • **In-Repo Rational Certificate (`certified_linear_strict_gap_sigma_quarter`):**
    Certifies the explicit rational numerical bounds at the benchmark point
    `σ = 1/4`, `δ = 1/20` (matching `validate/edge_w1_route_f_closure.py` and
    `validate/edge_w2_w3_closure.py`).
-/
import Mathlib
import SingularBlowup.PointSpectrum
import SingularBlowup.SlitSchurMargin

namespace SingularBlowup.LinearStrictGapClosure

open SingularBlowup.PointSpectrum
open SingularBlowup.SlitSchurMargin

noncomputable section

/-! ## 1. W1 Synthesis: Uniform Vertical-Line Resolvent Bound -/

/-- Total Route-F resolvent operator norm upper envelope at vertical frequency `ω = Im λ`
and spectral-gap margin `δ = Re λ + 2σ > 0`:
  `M_tot(ω) = (M_lip * C_ker + C_far) / δ + C_def / (1 + |ω|)`. -/
def routeFResolventEnvelope (M_lip C_ker C_far C_def delta omega : ℝ) : ℝ :=
  (M_lip * C_ker + C_far) / delta + C_def / (1 + |omega|)

/-- **W1 Master Theorem (Uniform Vertical-Line Resolvent Bound).**
Given the diagonal modulus-commutator bound `diagTerm ≤ M_lip * C_ker / δ` (from
`modulus_commutator_cancels_cauchy_pole`), the far-band Schur bound `farTerm ≤ C_far / δ`
(from `lem_stir_schur_four_quadrant_margin`), and the rank-one defect ω-cancellation
`defectTerm ≤ C_def / (1 + |ω|)` (from `slit_defect_rank_one_margin`), the total
resolvent norm at `λ = (-2σ + δ) + iω` is bounded uniformly for all `ω ∈ ℝ` by
`(M_lip * C_ker + C_far) / δ + C_def`. -/
theorem w1_uniform_resolvent_bound
    (M_lip C_ker C_far C_def delta omega diagTerm farTerm defectTerm : ℝ)
    (hCdef : 0 ≤ C_def)
    (h_diag : diagTerm ≤ M_lip * C_ker / delta)
    (h_far : farTerm ≤ C_far / delta)
    (h_def : defectTerm ≤ C_def / (1 + |omega|)) :
    diagTerm + farTerm + defectTerm ≤
      routeFResolventEnvelope M_lip C_ker C_far C_def delta omega ∧
    routeFResolventEnvelope M_lip C_ker C_far C_def delta omega ≤
      (M_lip * C_ker + C_far) / delta + C_def := by
  have h_env :
      diagTerm + farTerm + defectTerm ≤
        routeFResolventEnvelope M_lip C_ker C_far C_def delta omega := by
    unfold routeFResolventEnvelope
    have h_sum_div :
        M_lip * C_ker / delta + C_far / delta = (M_lip * C_ker + C_far) / delta := by
      ring
    linarith
  have h_unif :
      routeFResolventEnvelope M_lip C_ker C_far C_def delta omega ≤
        (M_lip * C_ker + C_far) / delta + C_def := by
    unfold routeFResolventEnvelope
    have h_denom : (1 : ℝ) ≤ 1 + |omega| := by
      have : (0 : ℝ) ≤ |omega| := abs_nonneg omega
      linarith
    have h_div : C_def / (1 + |omega|) ≤ C_def / 1 :=
      div_le_div_of_nonneg_left hCdef zero_lt_one h_denom
    linarith
  exact ⟨h_env, h_unif⟩

/-! ## 2. W2 Synthesis: Left-Strip Mixed Laplace Bridge + LEM-08 Uniqueness -/

/-- **W2 Master Theorem (Left-Strip Analyticity & Laurent Uniqueness without Edge Bootstrap).**
For any convective weight `σ ∈ (0, 1/2)`, spectral parameter with `δ = Re λ + 2σ > 0`,
and strip coordinate `σ' ∈ [σ - 1/2, σ]`:
1. The edge half-line `z ≤ 0` embeds `L²_σ((-∞, 0]) ⊆ L²_{σ'}((-∞, 0])` unconditionally;
2. The far-field half-line `z ≥ 0` embeds `KG ∈ L²_{σ-1/2}([0, ∞)) ⊆ L²_{σ'}([0, ∞))`;
3. All three Gamma factors in `(λ + 2s) X(s)` have strictly positive real parts (`≥ δ`, `> 1`, `> 1/2`);
4. Any integer Laurent mode support `[kMin, kMax]` compatible with the Stirling `L²` budget
   `3/4` collapses to `(0, 0)` and cannot kill the `s = 0` edge pole. -/
theorem w2_left_strip_bridge_and_uniqueness
    (sigma sigmap re_lam delta z_edge z_far : ℝ)
    (kMax kMin : ℤ)
    (h_sigma_pos : 0 < sigma)
    (h_sigma_half : sigma < 1 / 2)
    (h_delta_def : delta = re_lam + 2 * sigma)
    (h_delta_pos : 0 < delta)
    (h_lo : sigma - 1 / 2 ≤ sigmap)
    (h_hi : sigmap ≤ sigma)
    (hz_edge : z_edge ≤ 0)
    (hz_far : 0 ≤ z_far)
    (hkMax : (kMax : ℝ) < l2Budget)
    (hkMin : |(kMin : ℝ)| < l2Budget)
    (hkOrd : kMin ≤ kMax) :
    (-2 * sigmap * z_edge ≤ -2 * sigma * z_edge) ∧
    (-2 * sigmap * z_far ≤ -2 * (sigma - 1 / 2) * z_far) ∧
    (0 < 2 * sigmap + re_lam + 1 ∧ 1 < 2 * sigmap + 2 ∧ 1 / 2 < 1 - sigmap) ∧
    (kMax = 0 ∧ kMin = 0 ∧ ¬ ((1 : ℤ) ≤ laurentDegree kMax kMin)) := by
  have h_edge := left_strip_edge_exponent_le sigma sigmap z_edge hz_edge h_hi
  have h_far := left_strip_far_exponent_le sigma sigmap z_far hz_far h_lo
  have h_floors :=
    left_strip_gamma_real_floors sigma sigmap re_lam delta
      h_sigma_pos h_sigma_half h_delta_def h_lo h_hi
  have h_pos1 : 0 < 2 * sigmap + re_lam + 1 := lt_of_lt_of_le h_delta_pos h_floors.1
  have h_modes := only_constant_mode_l2_admissible kMax kMin hkMax hkMin hkOrd
  have h_nokill := no_l2_and_edge_kill kMax kMin hkMax hkMin hkOrd
  exact ⟨h_edge, h_far, ⟨h_pos1, h_floors.2.1, h_floors.2.2⟩, ⟨h_modes.1, h_modes.2, h_nokill⟩⟩

/-! ## 3. W3 arithmetic: the frequency pincer (its semigroup application is WITHDRAWN, 2026-11-21) -/

/-- **W3 High-Frequency Denominator Dichotomy.**
From `transport_coupling_frequency_pincer` (`|ω| ≤ |2τ + ω| + |2τ|`), at every
frequency `τ ∈ ℝ` and spectral height `ω = Im λ`, the sum of the regular-shifted
transport and coupling denominators `1 + |2τ + ω| + |2τ|` is at least `1 + |ω|`,
and the maximum of `|2τ + ω|` and `|2τ|` is at least `|ω| / 2`. -/
theorem w3_frequency_pincer_dichotomy (tau omega : ℝ) :
    |omega| / 2 ≤ max |2 * tau + omega| |2 * tau| ∧
    1 + |omega| ≤ 1 + |2 * tau + omega| + |2 * tau| := by
  have h_pincer := transport_coupling_frequency_pincer tau omega
  have h_le_left : |2 * tau + omega| ≤ max |2 * tau + omega| |2 * tau| := le_max_left _ _
  have h_le_right : |2 * tau| ≤ max |2 * tau + omega| |2 * tau| := le_max_right _ _
  constructor
  · linarith
  · linarith

/-- **W3 arithmetic (pointwise amplitude bound).**
Because `1 + |ω| ≤ 1 + |2τ + ω| + |2τ|` holds for all `(τ, ω) ∈ ℝ²`, any real number
`pertAmp` bounded by `C_pert / (1 + |2τ + ω| + |2τ|)` is bounded by `C_pert / (1 + |ω|)`.

This is all that is proved. It was intended to give `O((1 + |ω|)^{-1})` decay of
`‖R(λ) - R₀(λ)‖` along `Re λ = -2σ + δ`; that application is WITHDRAWN (2026-11-21), because
the hypothesis `h_amp` is not established as an operator bound on `L²_σ`
(see `docs/desk/w3_gearhart_pruss.md`, Correction). -/
theorem w3_perturbation_vertical_decay
    (C_pert tau omega pertAmp : ℝ)
    (hC : 0 ≤ C_pert)
    (h_amp : pertAmp ≤ C_pert / (1 + |2 * tau + omega| + |2 * tau|)) :
    pertAmp ≤ C_pert / (1 + |omega|) := by
  have h_pincer := (w3_frequency_pincer_dichotomy tau omega).2
  have h_pos : (0 : ℝ) < 1 + |omega| := by
    have : (0 : ℝ) ≤ |omega| := abs_nonneg omega
    linarith
  have h_div :
      C_pert / (1 + |2 * tau + omega| + |2 * tau|) ≤ C_pert / (1 + |omega|) :=
    div_le_div_of_nonneg_left hC h_pos h_pincer
  exact le_trans h_amp h_div

/-! ## 4. Rational arithmetic at the sample point `σ = 1/4`, `δ = 1/20` (not a certificate) -/

/-- Structure packaging constants for the arithmetic below at a point `(σ, δ)`.
The name is historical: nothing is certified about the operator by building a term of it. -/
structure CertifiedLinearGapWitness where
  sigma : ℝ
  delta : ℝ
  M_lip : ℝ
  C_ker : ℝ
  C_far : ℝ
  C_def : ℝ
  M_sup : ℝ
  sigma_pos : 0 < sigma
  sigma_lt_half : sigma < 1 / 2
  delta_pos : 0 < delta
  delta_lt_two_sigma : delta < 2 * sigma
  M_lip_nonneg : 0 ≤ M_lip
  C_ker_nonneg : 0 ≤ C_ker
  C_far_nonneg : 0 ≤ C_far
  C_def_nonneg : 0 ≤ C_def
  envelope_le_M_sup : (M_lip * C_ker + C_far) / delta + C_def ≤ M_sup

/-- **Sample constants at `σ = 1/4`, `δ = 1/20` (name kept; not a certified spectral gap).**

Caution (2026-11-21). These constants are illustrative and are NOT the constants of the
manuscript's proof: there the kernel constant is `C_ker = 4π/(1 - e^{-4π}) ≈ 12.57`, and the
diagonal bound of Lemma A.2 carries the factor `L e^{L}` with `L = L(1/4, 1/20) ≈ 54.76`, far
above the `1100` used in `h_diag` below. The "observed" figures are samples from the gate; a
wider scan finds `sup |∂_τ log|X|| ≈ 26.4` (near `τ ≈ -ω/2`), not `16.2276` — still below `55`.
So `M_sup = 1123` is an arithmetic consequence of the hypotheses below and of nothing else.

Original description (2026-11-20): matches the analytic Digamma bound and numerical evaluation
of `validate/edge_w1_route_f_closure.py` at `σ = 1/4`, `δ = 1/20` (`Re λ = -0.45`):
  • `σ = 1/4` (convective gap `2σ = 1/2`), `δ = 1/20` (certified half-plane `Re λ ≥ -9/20`),
  • `M_lip = 55` (analytic ceiling `L(1/4, 1/20) = 3π + 40 + 4 + 4/3 ≈ 54.758 < 55`; observed max `16.2276`),
  • `C_ker = 1` (split Cauchy diagonal kernel constant; observed commutator max `21.7818 ≤ 1100`),
  • `C_far = 1` (off-diagonal Schur bound `1 / (1/20) = 20`; observed max `0.34337`),
  • `C_def = 3` (rank-one defect sup `‖X‖₂ ‖K_slit‖₂ = 2.0689 < 3`, decaying to `0.0096` at `ω = 50`),
  • `M_sup = 1123` (`(55 * 1 + 1) / (1/20) + 3 = 1120 + 3 = 1123`). -/
def certifiedWitnessSigmaQuarter : CertifiedLinearGapWitness where
  sigma := 1 / 4
  delta := 1 / 20
  M_lip := 55
  C_ker := 1
  C_far := 1
  C_def := 3
  M_sup := 1123
  sigma_pos := by norm_num
  sigma_lt_half := by norm_num
  delta_pos := by norm_num
  delta_lt_two_sigma := by norm_num
  M_lip_nonneg := by norm_num
  C_ker_nonneg := by norm_num
  C_far_nonneg := by norm_num
  C_def_nonneg := by norm_num
  envelope_le_M_sup := by norm_num

/-- **Arithmetic at `σ = 1/4`, `δ = 1/20` (name kept; this is not a certified spectral gap).**
For the constants of `certifiedWitnessSigmaQuarter` and every real `ω`: any three real numbers
`diagTerm`, `farTerm`, `defectTerm` satisfying the three hypothesised bounds sum to at most
`M_sup = 1123`, and `-2σ + δ = -9/20`.

The hypotheses `h_diag`, `h_far`, `h_def` stand for operator-norm bounds on the line
`Re λ = -9/20` that are argued on paper in the manuscript (Appendix A) and are NOT proved here.
Nothing about the resolvent, the spectrum or a semigroup follows from this theorem alone. -/
theorem certified_linear_strict_gap_sigma_quarter
    (omega diagTerm farTerm defectTerm : ℝ)
    (h_diag : diagTerm ≤ certifiedWitnessSigmaQuarter.M_lip *
      certifiedWitnessSigmaQuarter.C_ker / certifiedWitnessSigmaQuarter.delta)
    (h_far : farTerm ≤ certifiedWitnessSigmaQuarter.C_far /
      certifiedWitnessSigmaQuarter.delta)
    (h_def : defectTerm ≤ certifiedWitnessSigmaQuarter.C_def / (1 + |omega|)) :
    -2 * certifiedWitnessSigmaQuarter.sigma + certifiedWitnessSigmaQuarter.delta = -9 / 20 ∧
    diagTerm + farTerm + defectTerm ≤ certifiedWitnessSigmaQuarter.M_sup := by
  have h_w1 :=
    w1_uniform_resolvent_bound
      certifiedWitnessSigmaQuarter.M_lip
      certifiedWitnessSigmaQuarter.C_ker
      certifiedWitnessSigmaQuarter.C_far
      certifiedWitnessSigmaQuarter.C_def
      certifiedWitnessSigmaQuarter.delta
      omega diagTerm farTerm defectTerm
      certifiedWitnessSigmaQuarter.C_def_nonneg
      h_diag h_far h_def
  have h_gap :
      -2 * certifiedWitnessSigmaQuarter.sigma + certifiedWitnessSigmaQuarter.delta = -9 / 20 := by
    simp [certifiedWitnessSigmaQuarter]
    norm_num
  have h_tot : diagTerm + farTerm + defectTerm ≤ certifiedWitnessSigmaQuarter.M_sup :=
    le_trans (le_trans h_w1.1 h_w1.2) certifiedWitnessSigmaQuarter.envelope_le_M_sup
  exact ⟨h_gap, h_tot⟩

end
end SingularBlowup.LinearStrictGapClosure
