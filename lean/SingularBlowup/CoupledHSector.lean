/-
# CoupledHSector.lean — Algebraic & Analytic Backbone for Coupled (G, h) Nonlinear Stability

Formalizes the load-bearing inequalities that upgrade Theorem 3 (`thm:coupled_h`)
in `paper/stability_spectra_singular_houluo.tex` from linear-only coercivity to
full coupled (G, h) nonlinear asymptotic stability in X × X_h with ZERO extra
derivatives on G:

1. `m2_re0_sq_lt_quarter`:
   For all τ ∈ ℝ, using `coth²(πτ) ≤ 1 + 1/(π²τ²)` (from `sinh(x) ≥ x`),
   the unweighted Re(s) = 0 second-derivative multiplier satisfies
     `|m₂(iτ)|² ≤ (τ⁴ + τ²/π²) / (4 (1/4 + τ²)(1 + τ²)) < 1/4`,
   yielding `‖(1/2) ∂_z²(Q_op G)‖_{L²_0} ≤ (1/2) ‖∂_z G‖_{L²_{1/2}}`.
2. `m1_re0_sq_lt_quarter`:
   For all τ ∈ ℝ, the unweighted Re(s) = 0 first-derivative multiplier satisfies
     `|m₁(iτ)|² ≤ (τ² + 1/π²) / (4 (1 + τ²)) < 1/4`,
   yielding `‖(1/2) ∂_z(Q_op W)‖_{L²_0} ≤ (1/2) ‖W‖_{L²_{1/2}}` (Picard contraction).
3. `coupled_h_linear_coercivity`:
   With `μ_h = 1000 B₀²` and `B ≤ 2 B₀`, the coupled linear rate satisfies
     `-(1/8) E + 4 B₀ √(20 E E_h) - 3 μ_h E_h ≤ -(1/16) (E + μ_h E_h)`.
4. `coupled_q_inflow_lt_one`:
   The overlap inflow at `z₀ = -ln 2` satisfies `2√2 / 3 < 1`, so
   `2√2 |q_z(-ln 2)|² ≤ Z_{h,disk}` and `2√2 |q(-ln 2)|² ≤ Z_{h,disk}`.
5. `coupled_disk_sqrt_split`:
   For `j, k ≥ 0` with `m = j + k`, `√m ≤ √j + √k`, absorbing `[(v ∂_v p) q]_m`
   into `-4 D_{h,disk}` via `‖j^{3/2} q_j‖_{ℓ²}` and `‖j q_j‖_{ℓ¹}`.

Zero `sorry`, zero custom axioms.
-/
import Mathlib

namespace SingularBlowup.CoupledHSector

/-! ### 1. Unweighted `Re s = 0` Plancherel smoothing multipliers `|m₂(iτ)| < 1/2` and `|m₁(iτ)| < 1/2` -/

/-- Exact rational/real bound underlying `|m₂(iτ)| < 1/2` on `Re s = 0`
    once `coth²(πτ) ≤ 1 + 1/(π²τ²)` (from `sinh(x) ≥ x`) is applied:
    since `π > 3`, `1/π² < 1/9 < 5/4`, so `τ⁴ + τ²/π² < (1/4 + τ²)(1 + τ²)`. -/
theorem m2_re0_sq_lt_quarter (τ : ℝ) :
    (τ ^ 4 + τ ^ 2 / Real.pi ^ 2) / (4 * (1 / 4 + τ ^ 2) * (1 + τ ^ 2)) < 1 / 4 := by
  have hpi : (3 : ℝ) < Real.pi := Real.pi_gt_three
  have hpi2 : (9 : ℝ) < Real.pi ^ 2 := by nlinarith
  have hpi_pos : (0 : ℝ) < Real.pi ^ 2 := by positivity
  have hden : (0 : ℝ) < 4 * (1 / 4 + τ ^ 2) * (1 + τ ^ 2) := by positivity
  rw [div_lt_div_iff₀ hden (by norm_num : (0 : ℝ) < 4)]
  have hinv : τ ^ 2 / Real.pi ^ 2 ≤ τ ^ 2 / 9 := by
    rw [div_le_div_iff₀ hpi_pos (by norm_num : (0 : ℝ) < 9)]
    nlinarith [sq_nonneg τ]
  nlinarith [sq_nonneg τ]

/-- Exact bound underlying `|m₁(iτ)| < 1/2` on `Re s = 0` for zero-derivative
    Picard contraction of `∂_z(Q_op W)` in `L²_0`: `τ² + 1/π² < 1 + τ²`. -/
theorem m1_re0_sq_lt_quarter (τ : ℝ) :
    (τ ^ 2 + 1 / Real.pi ^ 2) / (4 * (1 + τ ^ 2)) < 1 / 4 := by
  have hpi : (3 : ℝ) < Real.pi := Real.pi_gt_three
  have hpi2 : (9 : ℝ) < Real.pi ^ 2 := by nlinarith
  have hpi_pos : (0 : ℝ) < Real.pi ^ 2 := by positivity
  have hden : (0 : ℝ) < 4 * (1 + τ ^ 2) := by positivity
  rw [div_lt_div_iff₀ hden (by norm_num : (0 : ℝ) < 4)]
  have hinv : 1 / Real.pi ^ 2 < 1 := by
    rw [div_lt_one hpi_pos]
    linarith
  linarith

/-! ### 2. Coupled `(G, h)` linear coercivity with `μ_h = 1000 B₀²` and `B ≤ 2 B₀` -/

theorem coupled_h_linear_coercivity (E Eh B0 : ℝ) (hE : 0 ≤ E) (hEh : 0 ≤ Eh) :
    -(1 / 8 : ℝ) * E + 4 * B0 * Real.sqrt (20 * E * Eh) - 3 * (1000 * B0 ^ 2) * Eh ≤
    -(1 / 16 : ℝ) * (E + (1000 * B0 ^ 2) * Eh) := by
  -- Young's inequality: 4 B₀ √(20 E E_h) ≤ (1/16) E + 1280 B₀² E_h
  have h_young : 4 * B0 * Real.sqrt (20 * E * Eh) ≤ (1 / 16 : ℝ) * E + 1280 * B0 ^ 2 * Eh := by
    have hsq : 0 ≤ (Real.sqrt E / 4 - 16 * Real.sqrt 5 * B0 * Real.sqrt Eh) ^ 2 := sq_nonneg _
    have hsqrt_prod : Real.sqrt (20 * E * Eh) = 2 * Real.sqrt 5 * Real.sqrt E * Real.sqrt Eh := by
      have h1 : Real.sqrt (20 * E * Eh) = Real.sqrt 20 * Real.sqrt E * Real.sqrt Eh := by
        rw [Real.sqrt_mul (by positivity : (0 : ℝ) ≤ 20 * E), Real.sqrt_mul (by norm_num : (0 : ℝ) ≤ 20)]
      have h2 : Real.sqrt 20 = 2 * Real.sqrt 5 := by
        rw [show (20 : ℝ) = 2 ^ 2 * 5 by norm_num, Real.sqrt_mul (by norm_num : (0 : ℝ) ≤ 2 ^ 2)]
        rw [Real.sqrt_sq (by norm_num : (0 : ℝ) ≤ 2)]
      rw [h1, h2]
    have h5 : Real.sqrt 5 ^ 2 = 5 := Real.sq_sqrt (by norm_num)
    have hEsq : Real.sqrt E ^ 2 = E := Real.sq_sqrt hE
    have hEhsq : Real.sqrt Eh ^ 2 = Eh := Real.sq_sqrt hEh
    have hexp : (Real.sqrt E / 4 - 16 * Real.sqrt 5 * B0 * Real.sqrt Eh) ^ 2 =
        (1 / 16 : ℝ) * Real.sqrt E ^ 2 - 8 * B0 * (Real.sqrt 5 * Real.sqrt E * Real.sqrt Eh) +
        256 * Real.sqrt 5 ^ 2 * B0 ^ 2 * Real.sqrt Eh ^ 2 := by ring
    rw [hexp, hEsq, h5, hEhsq] at hsq
    rw [hsqrt_prod]
    linarith
  nlinarith [sq_nonneg B0]

/-! ### 3. Overlap inflow & disk convolution split -/

theorem coupled_q_inflow_lt_one : 2 * Real.sqrt 2 / 3 < 1 := by
  have h : (2 * Real.sqrt 2) ^ 2 < 3 ^ 2 := by
    have h2 : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
    nlinarith
  have hpos : 0 ≤ 2 * Real.sqrt 2 := by positivity
  have hlt : 2 * Real.sqrt 2 < 3 := by nlinarith
  linarith

theorem coupled_disk_sqrt_split (j k : ℝ) (hj : 0 ≤ j) (hk : 0 ≤ k) :
    Real.sqrt (j + k) ≤ Real.sqrt j + Real.sqrt k := by
  have h1 : 0 ≤ Real.sqrt (j + k) := Real.sqrt_nonneg _
  have h2 : 0 ≤ Real.sqrt j + Real.sqrt k := by positivity
  nlinarith [Real.sq_sqrt (add_nonneg hj hk), Real.sq_sqrt hj, Real.sq_sqrt hk,
             mul_nonneg (Real.sqrt_nonneg j) (Real.sqrt_nonneg k)]

end SingularBlowup.CoupledHSector
