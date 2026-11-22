/-
# WeakEdgeManifolds.lean — algebraic & pointwise analytic backbone of the
  2026-11-21 global nonlinear weak-edge compatibility manifolds (Tracks A, B, C)

Sources:
  - `docs/findings_session_2026_11_21_nonlinear_cusp_manifold.md` (Track A, `NM-1..3`)
  - `docs/findings_session_2026_11_21_pure_step_nonlinear_manifold.md` (Track B, `PM-1..3`)
  - `docs/findings_session_2026_11_21_finite_sobolev_pure_step.md` (Track C, `FS-1..3`)

Machine-checked (zero-sorry) theorems in this file:
  1. `phase_cancellation_norm` — exact pointwise phase-cancellation identity for ANY nonzero
     complex multiplier `m₁, m₂ ∈ ℂ`:
       `‖I * (ξ - η) * (m₂ / m₁) - I * ξ * star (m₁ / m₂)‖ = |(ξ - η) * (‖m₂‖ / ‖m₁‖) - ξ * (‖m₁‖ / ‖m₂‖)|`.
     This eliminates the Stirling log-phase `Φ(ξ) ~ 2ξ ln(2|ξ|)` identically from the
     symmetric part of `A_b = 𝒩_g⁻¹ (b ∂_z) 𝒩_g` on every `H^k_σ`.
  2. `hasDerivAt_scalarEdgeD` — the closed-form weak-edge cusp amplitude
     `D(s) = -a₀ / (1 + a₀ s / 2)` solves `D'(s) = (1/2) D(s)²` with `D(0) = -a₀`.
  3. `strip_re_u_ge`, `strip_norm_exp_le_one`, `strip_norm_us_exp_le` — sharp uniform bounds
     on `u = exp(-(x + i y)/2)` for `|y| ≤ ρ < π` and `s ≥ 0`:
       `Re(u) ≥ e^{-x/2} cos(ρ/2) > 0`, `‖exp(-u s)‖ ≤ 1`, `‖u s exp(-u s)‖ ≤ 1 / (e cos(ρ/2))`.
  4. `hierarchical_u_dissipation_H3` & `hierarchical_u_dissipation_H4` — pointwise quadratic-form
     coercivity for `∂_z^m(u v)` (`u = e^{-z/2}`, `∂_z u = -u/2`) with weights `c_m = 4^{3-m}`
     on `H³` and `c_m = 6^{4-m}` on `H⁴`:
       `∑_{m=0}^k c_m u ((∂_z - 1/2)^m v) (∂_z^m v) ≥ (1/2) ∑_{m=0}^k c_m u (∂_z^m v)²`.
  5. `nishida_contraction_bound` — an ARITHMETIC fact only: with `γ = 1 / (20 C_T M_G)` and
     `a = 40 ε / γ`, a three-term expression is `< 1/2`. It does not prove any estimate, and its
     parametrisation (cone speed `a ∝ ε`, `2√2`) differs from Lemma A3 of
     `docs/theorem_track_a_cusp_weak_edge.md` (`a ∈ [ρ₀/12, ρ₀/6]`, `γ = 1/(40 C_rem M_G)`, `4√2`).
     Lemma A3's contraction has a gap (2026-11-22 post-f5106e9 review): its norm does not close.

SCOPE. Operator-level Plancherel/Mellin bounds, the Nishida iteration in Banach scales, and
distributional PDE equivalence across `X = 1` remain on paper (`[D]`) and in the Python
gates (`NM-1..3`, `PM-1..3`, `FS-1..3`).
-/
import Mathlib

namespace SingularBlowup.WeakEdgeManifolds

open Complex Real

noncomputable section

/-! ### 1. Exact Phase-Cancellation Identity (Track C, `FS-1`) -/

/-- **Exact Phase-Cancellation Identity.** For any nonzero complex numbers `m₁, m₂` and real
frequencies `ξ, η`, the oscillatory phases of `m₂ / m₁` and `conj(m₁ / m₂)` are identical and
factor out with modulus `1`, leaving a purely amplitude-ratio expression. -/
theorem phase_cancellation_norm (m₁ m₂ : ℂ) (hm₁ : m₁ ≠ 0) (hm₂ : m₂ ≠ 0) (ξ η : ℝ) :
    ‖Complex.I * ((ξ - η : ℝ) : ℂ) * (m₂ / m₁) - Complex.I * (ξ : ℂ) * star (m₁ / m₂)‖ =
      |(ξ - η) * (‖m₂‖ / ‖m₁‖) - ξ * (‖m₁‖ / ‖m₂‖)| := by
  have h1 : (‖m₁‖ : ℝ) ≠ 0 := norm_ne_zero_iff.mpr hm₁
  have h2 : (‖m₂‖ : ℝ) ≠ 0 := norm_ne_zero_iff.mpr hm₂
  have h1c : ((‖m₁‖ : ℝ) : ℂ) ≠ 0 := Complex.ofReal_ne_zero.mpr h1
  have h2c : ((‖m₂‖ : ℝ) : ℂ) ≠ 0 := Complex.ofReal_ne_zero.mpr h2
  -- Unit phase factor u = (m₂ / ‖m₂‖) * (conj m₁ / ‖m₁‖)
  set u : ℂ := (m₂ / (‖m₂‖ : ℂ)) * (star m₁ / (‖m₁‖ : ℂ)) with hu_def
  have hu_norm : ‖u‖ = 1 := by
    rw [hu_def, norm_mul, norm_div, norm_div, Complex.norm_real, Complex.norm_real,
        norm_star, Real.norm_of_nonneg (norm_nonneg m₂), Real.norm_of_nonneg (norm_nonneg m₁),
        div_self h2, div_self h1, mul_one]
  -- Express m₂ / m₁ = u * (‖m₂‖ / ‖m₁‖)
  have hm_mul_star₁ : m₁ * star m₁ = ((‖m₁‖ : ℝ) : ℂ) ^ 2 := by
    change m₁ * starRingEnd ℂ m₁ = _
    rw [Complex.mul_conj, Complex.normSq_eq_norm_sq, Complex.ofReal_pow]
  have hm_mul_star₂ : m₂ * star m₂ = ((‖m₂‖ : ℝ) : ℂ) ^ 2 := by
    change m₂ * starRingEnd ℂ m₂ = _
    rw [Complex.mul_conj, Complex.normSq_eq_norm_sq, Complex.ofReal_pow]
  have h_ratio₁ : m₂ / m₁ = u * (((‖m₂‖ / ‖m₁‖ : ℝ) : ℂ)) := by
    rw [hu_def, Complex.ofReal_div]
    have hmul : m₂ / m₁ = (m₂ * star m₁) / (m₁ * star m₁) := by
      rw [mul_div_mul_right _ _ (star_ne_zero.mpr hm₁)]
    rw [hmul, hm_mul_star₁]
    field_simp [h1c, h2c]
  -- Express star (m₁ / m₂) = u * (‖m₁‖ / ‖m₂‖)
  have h_ratio₂ : star (m₁ / m₂) = u * (((‖m₁‖ / ‖m₂‖ : ℝ) : ℂ)) := by
    rw [hu_def, Complex.ofReal_div, star_div₀]
    have hmul : star m₁ / star m₂ = (m₂ * star m₁) / (m₂ * star m₂) := by
      rw [mul_div_mul_left _ _ hm₂]
    rw [hmul, hm_mul_star₂]
    field_simp [h1c, h2c]
  have h_factor :
      Complex.I * ((ξ - η : ℝ) : ℂ) * (m₂ / m₁) - Complex.I * (ξ : ℂ) * star (m₁ / m₂) =
        (Complex.I * u) * (((ξ - η) * (‖m₂‖ / ‖m₁‖) - ξ * (‖m₁‖ / ‖m₂‖) : ℝ) : ℂ) := by
    rw [h_ratio₁, h_ratio₂]
    push_cast
    ring
  rw [h_factor, norm_mul, norm_mul, Complex.norm_I, hu_norm, one_mul, one_mul,
      Complex.norm_real, Real.norm_eq_abs]

/-! ### 2. Exact Scalar Weak-Edge Compatibility Law (Track A, `NM-1`) -/

/-- Exact scalar edge cusp law `D(s) = -a₀ / (1 + a₀ s / 2)`. -/
def scalarEdgeD (a₀ s : ℝ) : ℝ := -a₀ / (1 + a₀ * s / 2)

theorem scalarEdgeD_zero (a₀ : ℝ) : scalarEdgeD a₀ 0 = -a₀ := by
  simp [scalarEdgeD]

/-- `D(s)` solves the exact nonlinear compatibility ODE `D'(s) = (1/2) D(s)²`
(i.e. `D' + (1/2) a D = 0` with `a = -D`). -/
theorem hasDerivAt_scalarEdgeD (a₀ s : ℝ) (hden : 1 + a₀ * s / 2 ≠ 0) :
    HasDerivAt (scalarEdgeD a₀) ((1 / 2 : ℝ) * (scalarEdgeD a₀ s) ^ 2) s := by
  have hlin : HasDerivAt (fun t : ℝ => 1 + a₀ * t / 2) (a₀ / 2) s := by
    have h1 : HasDerivAt (fun t : ℝ => a₀ * t) a₀ s := by
      simpa using (hasDerivAt_id s).const_mul a₀
    have h2 : HasDerivAt (fun t : ℝ => a₀ * t / 2) (a₀ / 2) s := h1.div_const 2
    simpa using h2.const_add (1 : ℝ)
  have hq := (hasDerivAt_const s (-a₀)).div hlin hden
  convert hq using 1
  unfold scalarEdgeD
  field_simp [hden]
  ring

/-! ### 3. Strip-Analytic Semigroup Kernel Bounds (Tracks A & B, `NM-2`) -/

/-- Complex `u(x + i y) = exp(-(x + i y)/2)`. -/
def stripU (x y : ℝ) : ℂ := Complex.exp (-((x : ℂ) + (y : ℂ) * Complex.I) / 2)

theorem stripU_re (x y : ℝ) : (stripU x y).re = Real.exp (-x / 2) * Real.cos (y / 2) := by
  unfold stripU
  rw [Complex.exp_re]
  have hre : (-((x : ℂ) + (y : ℂ) * Complex.I) / 2).re = -x / 2 := by
    simp [Complex.div_ofNat]
  have him : (-((x : ℂ) + (y : ℂ) * Complex.I) / 2).im = -(y / 2) := by
    simp [Complex.div_ofNat]; ring
  rw [hre, him, Real.cos_neg]

theorem stripU_norm (x y : ℝ) : ‖stripU x y‖ = Real.exp (-x / 2) := by
  unfold stripU
  rw [Complex.norm_exp]
  congr 1
  have harg : -((x : ℂ) + (y : ℂ) * Complex.I) / 2 = ((-x / 2 : ℝ) : ℂ) + ((-y / 2 : ℝ) : ℂ) * Complex.I := by
    push_cast; ring
  rw [harg]
  simp

theorem strip_re_u_ge (x y ρ : ℝ) (hρ : ρ < Real.pi) (hy : |y| ≤ ρ) :
    Real.exp (-x / 2) * Real.cos (ρ / 2) ≤ (stripU x y).re := by
  rw [stripU_re]
  apply mul_le_mul_of_nonneg_left _ (Real.exp_pos _).le
  have hy2 : |y / 2| ≤ ρ / 2 := by
    rw [abs_div, abs_of_pos (by norm_num : (0 : ℝ) < 2)]
    linarith
  have hρ_nn : 0 ≤ ρ / 2 := le_trans (abs_nonneg (y / 2)) hy2
  have hρ_lt : ρ / 2 < Real.pi / 2 := by linarith
  rw [← Real.cos_abs (y / 2)]
  exact Real.strictAntiOn_cos.antitoneOn
    ⟨abs_nonneg _, by linarith⟩
    ⟨hρ_nn, by linarith⟩
    hy2

theorem strip_norm_exp_le_one (x y ρ s : ℝ) (hρ : ρ < Real.pi) (hy : |y| ≤ ρ) (hs : 0 ≤ s) :
    ‖Complex.exp (-stripU x y * (s : ℂ))‖ ≤ 1 := by
  rw [Complex.norm_exp]
  have hre_cos : 0 < Real.cos (ρ / 2) := by
    apply Real.cos_pos_of_mem_Ioo
    have hρ_nn : 0 ≤ ρ := le_trans (abs_nonneg y) hy
    constructor <;> linarith
  have hre : 0 ≤ (stripU x y).re := by
    have h1 := strip_re_u_ge x y ρ hρ hy
    have h2 : 0 ≤ Real.exp (-x / 2) * Real.cos (ρ / 2) :=
      mul_nonneg (Real.exp_pos _).le hre_cos.le
    linarith
  have hre_prod : (-stripU x y * (s : ℂ)).re = -((stripU x y).re * s) := by
    simp [Complex.mul_re]
  rw [hre_prod]
  exact Real.exp_le_one_iff.mpr (by nlinarith)

/-- Uniform bound `‖u s exp(-u s)‖ ≤ 1 / (e cos(ρ/2))` on the strip `|y| ≤ ρ < π`, `s ≥ 0`. -/
theorem strip_norm_us_exp_le (x y ρ s : ℝ) (hρ : ρ < Real.pi) (hy : |y| ≤ ρ) (hs : 0 ≤ s) :
    ‖stripU x y * (s : ℂ) * Complex.exp (-stripU x y * (s : ℂ))‖ ≤
      1 / (Real.exp 1 * Real.cos (ρ / 2)) := by
  have hρ_nn : 0 ≤ ρ := le_trans (abs_nonneg y) hy
  have hcos : 0 < Real.cos (ρ / 2) :=
    Real.cos_pos_of_mem_Ioo ⟨by linarith, by linarith⟩
  have hexp_pos : 0 < Real.exp (-x / 2) := Real.exp_pos _
  have hre_ge := strip_re_u_ge x y ρ hρ hy
  have hre_prod : (-stripU x y * (s : ℂ)).re = -((stripU x y).re * s) := by
    simp [Complex.mul_re]
  rw [norm_mul, norm_mul, stripU_norm, Complex.norm_real, Real.norm_of_nonneg hs,
      Complex.norm_exp, hre_prod]
  -- Let t = Real.exp (-x / 2) * Real.cos (ρ / 2) * s ≥ 0.
  -- Since t ≤ (stripU x y).re * s, exp(-re*s) ≤ exp(-t), and t * exp(-t) ≤ exp(-1) = 1 / exp 1.
  set t : ℝ := Real.exp (-x / 2) * Real.cos (ρ / 2) * s with ht_def
  have ht_le : t ≤ (stripU x y).re * s := by
    rw [ht_def]
    exact mul_le_mul_of_nonneg_right hre_ge hs
  have hexp_le : Real.exp (-((stripU x y).re * s)) ≤ Real.exp (-t) :=
    Real.exp_le_exp.mpr (neg_le_neg ht_le)
  have hte : t * Real.exp (-t) ≤ Real.exp (-1) := by
    -- Equivalent to t ≤ exp(t - 1), i.e. 1 + (t - 1) ≤ exp(t - 1)
    have h1 : t ≤ Real.exp (t - 1) := by
      have := Real.add_one_le_exp (t - 1)
      linarith
    have h2 : Real.exp (t - 1) * Real.exp (-t) = Real.exp (-1) := by
      rw [← Real.exp_add]; congr 1; ring
    calc t * Real.exp (-t) ≤ Real.exp (t - 1) * Real.exp (-t) :=
          mul_le_mul_of_nonneg_right h1 (Real.exp_pos _).le
      _ = Real.exp (-1) := h2
  have hlhs_eq : Real.exp (-x / 2) * s * Real.exp (-t) = (t * Real.exp (-t)) / Real.cos (ρ / 2) := by
    rw [ht_def]
    field_simp [hcos.ne']
  calc Real.exp (-x / 2) * s * Real.exp (-((stripU x y).re * s))
      ≤ Real.exp (-x / 2) * s * Real.exp (-t) :=
        mul_le_mul_of_nonneg_left hexp_le (mul_nonneg hexp_pos.le hs)
    _ = (t * Real.exp (-t)) / Real.cos (ρ / 2) := hlhs_eq
    _ ≤ Real.exp (-1) / Real.cos (ρ / 2) :=
        div_le_div_of_nonneg_right hte hcos.le
    _ = 1 / (Real.exp 1 * Real.cos (ρ / 2)) := by
        rw [Real.exp_neg, inv_eq_one_div, div_div]

/-! ### 4. Hierarchical `u`-Dissipation Coercivity on `H³` (Track C, `FS-2`) -/

/-- **Pointwise `H³` hierarchical `u`-dissipation coercivity.**
Since `∂_z(u) = -(1/2) u`, `∂_z^m(u v) = u (∂_z - 1/2)^m v`. With weights
`(c₀, c₁, c₂, c₃) = (64, 16, 4, 1)`, the quadratic form is strictly coercive with margin `1/2`
for every `u ≥ 0` and every `(v₀, v₁, v₂, v₃) ∈ ℝ⁴`. -/
theorem hierarchical_u_dissipation_H3 (u v₀ v₁ v₂ v₃ : ℝ) (hu : 0 ≤ u) :
    64 * (u * v₀) * v₀
      + 16 * (u * (v₁ - (1 / 2 : ℝ) * v₀)) * v₁
      + 4 * (u * (v₂ - v₁ + (1 / 4 : ℝ) * v₀)) * v₂
      + 1 * (u * (v₃ - (3 / 2 : ℝ) * v₂ + (3 / 4 : ℝ) * v₁ - (1 / 8 : ℝ) * v₀)) * v₃
      ≥ (1 / 2 : ℝ) * (64 * u * v₀ ^ 2 + 16 * u * v₁ ^ 2 + 4 * u * v₂ ^ 2 + 1 * u * v₃ ^ 2) := by
  -- Factor out u ≥ 0; the remaining quadratic form in (v₀, v₁, v₂, v₃) is an explicit sum of squares
  have hq :
      0 ≤ 32 * v₀ ^ 2 - 8 * v₀ * v₁ + 8 * v₁ ^ 2 + v₀ * v₂ - 4 * v₁ * v₂ + 2 * v₂ ^ 2
          - (1 / 8 : ℝ) * v₀ * v₃ + (3 / 4 : ℝ) * v₁ * v₃ - (3 / 2 : ℝ) * v₂ * v₃
          + (1 / 2 : ℝ) * v₃ ^ 2 := by
    have hsos :
        32 * v₀ ^ 2 - 8 * v₀ * v₁ + 8 * v₁ ^ 2 + v₀ * v₂ - 4 * v₁ * v₂ + 2 * v₂ ^ 2
            - (1 / 8 : ℝ) * v₀ * v₃ + (3 / 4 : ℝ) * v₁ * v₃ - (3 / 2 : ℝ) * v₂ * v₃
            + (1 / 2 : ℝ) * v₃ ^ 2
          = (1 / 2 : ℝ) * (v₃ - (3 / 2 : ℝ) * v₂ + (3 / 4 : ℝ) * v₁ - (1 / 8 : ℝ) * v₀) ^ 2
            + (7 / 8 : ℝ) * (v₂ - (23 / 14 : ℝ) * v₁ + (13 / 28 : ℝ) * v₀) ^ 2
            + (75 / 14 : ℝ) * (v₁ - (46 / 75 : ℝ) * v₀) ^ 2
            + (17873 / 600 : ℝ) * v₀ ^ 2 := by ring
    rw [hsos]
    positivity
  have hmul := mul_nonneg hu hq
  linarith

/-- **Pointwise `H⁴` hierarchical `u`-dissipation coercivity (Track C, `TCB-4`).**
With weights `(c₀, c₁, c₂, c₃, c₄) = (1296, 216, 36, 6, 1)` (`c_m = 6^{4-m}`), the quadratic form
for `∂_z^m(u v) = u (∂_z - 1/2)^m v` is strictly coercive with margin `1/2` for every `u ≥ 0`
and every `(v₀, v₁, v₂, v₃, v₄) ∈ ℝ⁵`. -/
theorem hierarchical_u_dissipation_H4 (u v₀ v₁ v₂ v₃ v₄ : ℝ) (hu : 0 ≤ u) :
    1296 * (u * v₀) * v₀
      + 216 * (u * (v₁ - (1 / 2 : ℝ) * v₀)) * v₁
      + 36 * (u * (v₂ - v₁ + (1 / 4 : ℝ) * v₀)) * v₂
      + 6 * (u * (v₃ - (3 / 2 : ℝ) * v₂ + (3 / 4 : ℝ) * v₁ - (1 / 8 : ℝ) * v₀)) * v₃
      + 1 * (u * (v₄ - 2 * v₃ + (3 / 2 : ℝ) * v₂ - (1 / 2 : ℝ) * v₁ + (1 / 16 : ℝ) * v₀)) * v₄
      ≥ (1 / 2 : ℝ) * (1296 * u * v₀ ^ 2 + 216 * u * v₁ ^ 2 + 36 * u * v₂ ^ 2
          + 6 * u * v₃ ^ 2 + 1 * u * v₄ ^ 2) := by
  have hq :
      0 ≤ 648 * v₀ ^ 2 - 108 * v₀ * v₁ + 108 * v₁ ^ 2 + 9 * v₀ * v₂ - 36 * v₁ * v₂ + 18 * v₂ ^ 2
          - (3 / 4 : ℝ) * v₀ * v₃ + (9 / 2 : ℝ) * v₁ * v₃ - 9 * v₂ * v₃ + 3 * v₃ ^ 2
          + (1 / 16 : ℝ) * v₀ * v₄ - (1 / 2 : ℝ) * v₁ * v₄ + (3 / 2 : ℝ) * v₂ * v₄
          - 2 * v₃ * v₄ + (1 / 2 : ℝ) * v₄ ^ 2 := by
    have hsos :
        648 * v₀ ^ 2 - 108 * v₀ * v₁ + 108 * v₁ ^ 2 + 9 * v₀ * v₂ - 36 * v₁ * v₂ + 18 * v₂ ^ 2
            - (3 / 4 : ℝ) * v₀ * v₃ + (9 / 2 : ℝ) * v₁ * v₃ - 9 * v₂ * v₃ + 3 * v₃ ^ 2
            + (1 / 16 : ℝ) * v₀ * v₄ - (1 / 2 : ℝ) * v₁ * v₄ + (3 / 2 : ℝ) * v₂ * v₄
            - 2 * v₃ * v₄ + (1 / 2 : ℝ) * v₄ ^ 2
          = (1 / 2 : ℝ) * (v₄ - 2 * v₃ + (3 / 2 : ℝ) * v₂ - (1 / 2 : ℝ) * v₁ + (1 / 16 : ℝ) * v₀) ^ 2
            + (v₃ - 3 * v₂ + (7 / 4 : ℝ) * v₁ - (5 / 16 : ℝ) * v₀) ^ 2
            + (63 / 8 : ℝ) * (v₂ - (11 / 7 : ℝ) * v₁ + (25 / 56 : ℝ) * v₀) ^ 2
            + (9561 / 112 : ℝ) * (v₁ - (7155 / 12748 : ℝ) * v₀) ^ 2
            + (63172875 / 101984 : ℝ) * v₀ ^ 2 := by ring
    rw [hsos]
    positivity
  have hmul := mul_nonneg hu hq
  linarith

/-! ### 5. Nishida–Ovsyannikov Contraction Closure (Track A, §3) -/

/-- Arithmetic only: with `γ = 1 / (20 C_T M_G)` and `a = 40 ε / γ`, the expression
`2ε/(γa) + 2√2 γ C_T M_G + C_rem ε/a` is `< 1/2` (assuming `C_rem / (800 C_T M_G) ≤ 1/10`).
This does NOT verify the coupled Nishida–Ovsyannikov estimate of Track A Lemma A3, whose weighted
norm does not close (2026-11-22 post-f5106e9 review), and its parametrisation differs from that lemma's. -/
theorem nishida_contraction_bound (ε C_T M_G C_rem : ℝ)
    (hε : 0 < ε) (hCT : 0 < C_T) (hMG : 0 < M_G)
    (hrem : C_rem / (800 * C_T * M_G) ≤ (1 / 10 : ℝ)) :
    let γ : ℝ := 1 / (20 * C_T * M_G)
    let a : ℝ := 40 * ε / γ
    2 * ε / (γ * a) + 2 * Real.sqrt 2 * γ * C_T * M_G + C_rem * ε / a < (1 / 2 : ℝ) := by
  intro γ a
  have hpos : 0 < C_T * M_G := mul_pos hCT hMG
  have hγ : 0 < γ := by positivity
  have ha : a = 800 * C_T * M_G * ε := by
    simp only [a, γ]
    field_simp [hpos.ne', hε.ne']
    ring
  have h1 : 2 * ε / (γ * a) = 1 / 20 := by
    rw [ha]
    simp only [γ]
    field_simp [hpos.ne', hε.ne']
    ring
  have h2 : 2 * Real.sqrt 2 * γ * C_T * M_G = Real.sqrt 2 / 10 := by
    simp only [γ]
    field_simp [hpos.ne']
    ring
  have h3 : C_rem * ε / a = C_rem / (800 * C_T * M_G) := by
    rw [ha]
    field_simp [hpos.ne', hε.ne']
  have hsqrt : Real.sqrt 2 < (3 / 2 : ℝ) := by
    rw [Real.sqrt_lt' (by norm_num)]
    norm_num
  rw [h1, h2, h3]
  linarith

end

end SingularBlowup.WeakEdgeManifolds
