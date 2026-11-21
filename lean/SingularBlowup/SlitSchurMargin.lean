/-
# SlitSchurMargin.lean — W1 Route F algebraic & exponential-margin backbone

Formalizes the four load-bearing inequalities behind the Route-F uniform
L²_σ resolvent bound for `𝓕 ↦ 𝒢(s) = X(s) [D(e^{4πis}) - D(1)]`
(`docs/desk/w1_plemelj_stirling_writeup.md`, Route F, 2026-11-20):

1. **`lem_stir_schur_four_quadrant_margin` (Lem Stir-Schur):**
   For every `(τ, τ', ω) ∈ ℝ³` (where `ω = Im λ`), the net Stirling + split-Cauchy
   exponential decay exponent
     `E_ω(τ, τ') = π(|τ + ω/2| - |τ' + ω/2|) + 2π(|τ| - |τ'|) + 4π · splitWing(τ, τ')`
   satisfies `π |τ - τ'| ≤ E_ω(τ, τ')` across all four sign quadrants of `(τ, τ')`
   and uniformly in `ω = Im λ`.

2. **`slit_defect_rank_one_margin` (Lem D-slit-CS rank-one ω-cancellation):**
   For the rank-one defect operator `(T_defect 𝓕)(τ) = -X(σ+iτ; λ) D(1)[𝓕]`,
   the combined exponent
     `(π|τ + ω/2| + 2π|τ|) + (2π|τ'| - π|τ' + ω/2|)`
   dominates the separable decay bound `π|τ| + π|τ'|` uniformly in `ω = Im λ ∈ ℝ`.

3. **`shifted_weight_submultiplicative` (ω-uniform polynomial prefactor bound):**
   For every `(τ, τ', ω) ∈ ℝ³`,
     `1 + |τ + ω/2| ≤ (1 + |τ' + ω/2|) * (1 + |τ - τ'|)`,
   so the shifted Stirling polynomial ratios depend only on the translation `|τ - τ'|`,
   uniformly in `ω = Im λ`.

4. **`modulus_commutator_cancels_cauchy_pole` (Lem Diag-Riesz unitary phase reduction):**
   After factoring out the unitary phase `u(τ; λ) = X/|X|` (`|u| = 1`), the real modulus
   ratio satisfies `||X(τ)|/|X(τ')| - 1| ≤ M_lip * |τ - τ'|`, which cancels the Cauchy pole
   `C_ker / |τ - τ'|` on the diagonal band `0 < |τ - τ'| < 1`.
-/
import Mathlib

namespace SingularBlowup.SlitSchurMargin

noncomputable section

/-- Wing exponent of the split-regularized Cauchy kernel `ker(τ, τ')`:
plain kernel on `τ' ≥ 0` gives exponential decay `4π max(0, τ' - τ)`,
and the one-subtraction kernel on `τ' < 0` gives `4π max(0, τ - τ')`. -/
def splitWingExponent (tau taup : ℝ) : ℝ :=
  if 0 ≤ taup then max 0 (taup - tau) else max 0 (tau - taup)

/-- Net exponential decay exponent of the dressed split Cauchy kernel
`|X(σ+iτ; λ) / X(σ+iτ'; λ)| · ker(τ, τ')` at `ω = Im λ`. -/
def dressedDecayExponent (tau taup omega : ℝ) : ℝ :=
  Real.pi * (|tau + omega / 2| - |taup + omega / 2|) +
    2 * Real.pi * (|tau| - |taup|) +
    4 * Real.pi * splitWingExponent tau taup

/-- Auxiliary triangle bound: `-(|a - b|) ≤ |a| - |b|`. -/
theorem neg_abs_sub_le_abs_sub_abs (a b : ℝ) : -|a - b| ≤ |a| - |b| := by
  have h : |b| - |a| ≤ |b - a| := abs_sub_abs_le_abs_sub b a
  have h2 : |b - a| = |a - b| := abs_sub_comm b a
  linarith

/-- **Core 4-quadrant linear inequality (unit π-free form).**
For all `τ, τ', ω ∈ ℝ`,
`|τ - τ'| ≤ (|τ + ω/2| - |τ' + ω/2|) + 2(|τ| - |τ'|) + 4 · splitWingExponent τ τ'`. -/
theorem four_quadrant_margin_unit (tau taup omega : ℝ) :
    |tau - taup| ≤
      (|tau + omega / 2| - |taup + omega / 2|) +
        2 * (|tau| - |taup|) +
        4 * splitWingExponent tau taup := by
  have h_shift : -|tau - taup| ≤ |tau + omega / 2| - |taup + omega / 2| := by
    have h1 := neg_abs_sub_le_abs_sub_abs (tau + omega / 2) (taup + omega / 2)
    have h2 : (tau + omega / 2) - (taup + omega / 2) = tau - taup := by ring
    simpa [h2] using h1
  have h_base : -|tau - taup| ≤ |tau| - |taup| := neg_abs_sub_le_abs_sub_abs tau taup
  unfold splitWingExponent
  split_ifs with h_taup
  · by_cases h_ord : taup ≤ tau
    · have h_tau_nonneg : 0 ≤ tau := le_trans h_taup h_ord
      have h_abs_tau : |tau| = tau := abs_of_nonneg h_tau_nonneg
      have h_abs_taup : |taup| = taup := abs_of_nonneg h_taup
      have h_abs_diff : |tau - taup| = tau - taup := abs_of_nonneg (sub_nonneg.mpr h_ord)
      have h_max : max 0 (taup - tau) = 0 := max_eq_left (by linarith)
      linarith
    · have h_lt : tau < taup := lt_of_not_ge h_ord
      have h_abs_diff : |tau - taup| = taup - tau := by
        rw [abs_sub_comm]
        exact abs_of_nonneg (by linarith)
      have h_max : max 0 (taup - tau) = taup - tau := max_eq_right (by linarith)
      linarith
  · have h_taup_neg : taup < 0 := lt_of_not_ge h_taup
    by_cases h_ord : tau ≤ taup
    · have h_tau_neg : tau < 0 := lt_of_le_of_lt h_ord h_taup_neg
      have h_abs_tau : |tau| = -tau := abs_of_neg h_tau_neg
      have h_abs_taup : |taup| = -taup := abs_of_neg h_taup_neg
      have h_abs_diff : |tau - taup| = taup - tau := by
        rw [abs_sub_comm]
        exact abs_of_nonneg (by linarith)
      have h_max : max 0 (tau - taup) = 0 := max_eq_left (by linarith)
      linarith
    · have h_lt : taup < tau := lt_of_not_ge h_ord
      have h_abs_diff : |tau - taup| = tau - taup := abs_of_nonneg (by linarith)
      have h_max : max 0 (tau - taup) = tau - taup := max_eq_right (by linarith)
      linarith

/-- **Lem Stir-Schur (four-quadrant + all-ω exponential decay margin).**
For all `τ, τ', ω ∈ ℝ`, the dressed split Cauchy kernel decay exponent satisfies
`π |τ - τ'| ≤ dressedDecayExponent τ τ' ω`. -/
theorem lem_stir_schur_four_quadrant_margin (tau taup omega : ℝ) :
    Real.pi * |tau - taup| ≤ dressedDecayExponent tau taup omega := by
  unfold dressedDecayExponent
  have h := four_quadrant_margin_unit tau taup omega
  have hpi : 0 < Real.pi := Real.pi_pos
  have hmul := mul_le_mul_of_nonneg_left h (le_of_lt hpi)
  linarith

/-- **Lem D-slit-CS (rank-one defect ω-cancellation and separable decay margin).**
In `(T_defect 𝓕)(τ) = -X(σ+iτ; λ) D(1)[𝓕]`, the exponential decay exponent of `|X(τ)|` is
`π|τ + ω/2| + 2π|τ|` and that of `|K_slit(τ')|` is `2π|τ'| - π|τ' + ω/2|`.
Their sum dominates `π|τ| + π|τ'|` uniformly for all `τ, τ', ω ∈ ℝ`. -/
theorem slit_defect_rank_one_margin (tau taup omega : ℝ) :
    Real.pi * |tau| + Real.pi * |taup| ≤
      (Real.pi * |tau + omega / 2| + 2 * Real.pi * |tau|) +
        (2 * Real.pi * |taup| - Real.pi * |taup + omega / 2|) := by
  have h1 : |omega / 2| - |tau| ≤ |tau + omega / 2| := by
    have h_sub := abs_sub_abs_le_abs_sub (omega / 2) (-tau)
    have h_eq : (omega / 2) - (-tau) = tau + omega / 2 := by ring
    rw [h_eq, abs_neg] at h_sub
    exact h_sub
  have h2 : |taup + omega / 2| - |taup| ≤ |omega / 2| := by
    have h_sub := abs_sub_abs_le_abs_sub (taup + omega / 2) taup
    have h_eq : (taup + omega / 2) - taup = omega / 2 := by ring
    simpa [h_eq] using h_sub
  have h_unit : |tau| + |taup| ≤ (|tau + omega / 2| + 2 * |tau|) + (2 * |taup| - |taup + omega / 2|) := by
    linarith
  have hpi : 0 < Real.pi := Real.pi_pos
  have hmul := mul_le_mul_of_nonneg_left h_unit (le_of_lt hpi)
  linarith

/-- **ω-uniform submultiplicativity of shifted polynomial Stirling weights.**
For all `τ, τ', ω ∈ ℝ`, `1 + |τ + ω/2| ≤ (1 + |τ' + ω/2|) * (1 + |τ - τ'|)`. -/
theorem shifted_weight_submultiplicative (tau taup omega : ℝ) :
    1 + |tau + omega / 2| ≤ (1 + |taup + omega / 2|) * (1 + |tau - taup|) := by
  have h_tri : |tau + omega / 2| - |taup + omega / 2| ≤ |tau - taup| := by
    have h := abs_sub_abs_le_abs_sub (tau + omega / 2) (taup + omega / 2)
    have h_eq : (tau + omega / 2) - (taup + omega / 2) = tau - taup := by ring
    simpa [h_eq] using h
  have h_prod : 0 ≤ |taup + omega / 2| * |tau - taup| := mul_nonneg (abs_nonneg _) (abs_nonneg _)
  nlinarith

/-- **Lem Diag-Riesz (modulus commutator cancels the Cauchy pole).**
On the diagonal band `0 < d = |τ - τ'|`, a Lipschitz bound `|Δ| ≤ M_lip * d` on the
real modulus ratio `Δ = |X(τ)|/|X(τ')| - 1` cancels the Cauchy singularity `K ≤ C_ker / d`
and yields the uniform bound `|Δ| * K / δ ≤ M_lip * C_ker / δ`. -/
theorem modulus_commutator_cancels_cauchy_pole
    (M_lip C_ker delta d abs_Delta K : ℝ)
    (hd : 0 < d)
    (hdelta : 0 < delta)
    (hK_nonneg : 0 ≤ K)
    (hM_nonneg : 0 ≤ M_lip)
    (h_lip : abs_Delta ≤ M_lip * d)
    (h_ker : K ≤ C_ker / d) :
    abs_Delta * K / delta ≤ M_lip * C_ker / delta := by
  have h_mul1 : abs_Delta * K ≤ (M_lip * d) * (C_ker / d) :=
    mul_le_mul h_lip h_ker hK_nonneg (mul_nonneg hM_nonneg (le_of_lt hd))
  have h_cancel : (M_lip * d) * (C_ker / d) = M_lip * C_ker := by
    field_simp [ne_of_gt hd]
  rw [h_cancel] at h_mul1
  exact div_le_div_of_nonneg_right h_mul1 (le_of_lt hdelta)


/-- **W2 Left-Strip Paley–Wiener Split (edge half-line monotonicity).**
For `z ≤ 0` and `σ' ∈ [σ - 1/2, σ]`, `-2 * σ' * z ≤ -2 * σ * z`, so `e^{-2σ'z} ≤ e^{-2σz}`
and `L²_σ((-∞, 0]) ⊆ L²_{σ'}((-∞, 0])` unconditionally without an edge bootstrap. -/
theorem left_strip_edge_exponent_le (sigma sigmap z : ℝ)
    (hz : z ≤ 0) (h_hi : sigmap ≤ sigma) :
    -2 * sigmap * z ≤ -2 * sigma * z := by
  nlinarith

/-- **W2 Left-Strip Paley–Wiener Split (far-field half-line monotonicity).**
For `z ≥ 0` and `σ' ∈ [σ - 1/2, σ]`, `-2 * σ' * z ≤ -2 * (σ - 1/2) * z`, so
`KG = e^{-z/2} T_c G ∈ L²_{σ-1/2}([0, ∞))` embeds into `L²_{σ'}([0, ∞))` across the left strip. -/
theorem left_strip_far_exponent_le (sigma sigmap z : ℝ)
    (hz : 0 ≤ z) (h_lo : sigma - 1 / 2 ≤ sigmap) :
    -2 * sigmap * z ≤ -2 * (sigma - 1 / 2) * z := by
  nlinarith

/-- **W2 Isometric weight shift for `K = e^{-z/2} T_c`.**
The exponent `-z` from `|e^{-z/2}|²` and `-2(σ - 1/2)z` sum identically to `-2σz`. -/
theorem shift_half_weight_exponent_eq (sigma z : ℝ) :
    -z + (-2 * (sigma - 1 / 2) * z) = -2 * sigma * z := by
  ring

/-- **W2 Left-strip Gamma real-part positivity floors.**
Across the entire closed left strip `σ' ∈ [σ - 1/2, σ]` with `σ ∈ (0, 1/2)` and
`δ = Re λ + 2σ > 0`, the arguments of `Γ(2s + λ + 1)`, `Γ(2s + 2)`, and `Γ(1 - s)` in
`(λ + 2s) X(s) = Γ(2s + λ + 1) Γ(2s + 2) Γ(s) Γ(1 - s)` remain strictly inside the
pole-free right half-plane (`≥ δ`, `≥ 1`, and `≥ 1/2` respectively). -/
theorem left_strip_gamma_real_floors (sigma sigmap re_lam delta : ℝ)
    (h_sigma_pos : 0 < sigma)
    (h_sigma_half : sigma < 1 / 2)
    (h_delta_def : delta = re_lam + 2 * sigma)
    (h_lo : sigma - 1 / 2 ≤ sigmap)
    (h_hi : sigmap ≤ sigma) :
    delta ≤ 2 * sigmap + re_lam + 1 ∧
      1 < 2 * sigmap + 2 ∧
      1 / 2 < 1 - sigmap := by
  refine ⟨by linarith, by linarith, by linarith⟩

/-- **W3 High-frequency transport-vs-coupling separation pincer.**
For every frequency `τ` and vertical spectral parameter `ω = Im λ`,
`|ω| ≤ |2τ + ω| + |2τ|`, so at least one of the transport denominator `|λ + 2s| ≥ |2τ + ω|`
or the coupling denominator `|2s + 1| ≥ |2τ|` is `≥ |ω| / 2`. -/
theorem transport_coupling_frequency_pincer (tau omega : ℝ) :
    |omega| ≤ |2 * tau + omega| + |2 * tau| := by
  have h := abs_sub_abs_le_abs_sub omega (-(2 * tau))
  have h_eq : omega - -(2 * tau) = 2 * tau + omega := by ring
  rw [h_eq, abs_neg] at h
  linarith

end
end SingularBlowup.SlitSchurMargin
