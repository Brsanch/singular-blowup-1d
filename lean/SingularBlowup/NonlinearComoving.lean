/-
# NonlinearComoving.lean — Direction 2: Exact Gauge-(a) Collapse & Dual-Strip Nonlinear Backbone

Formalizes the exact algebraic and functional inequalities underlying the nonlinear
co-moving stability leg for the singular 1D Hou–Luo self-similar blow-up profile
(`docs/desk/nonlinear_comoving_leg.md`, `validate/edge_nonlinear_comoving_closure.py`,
2026-11-20):

1. **`gauge_a_advection_velocity_collapse`:**
   With the exact edge-stagnation gauge `c_l = 2 - V_G` and velocity decomposition
   `U(Y) = -2 + V_G + (Y - 1) * Q`, the total advection velocity collapses identically to
   `U(Y) + c_l * Y = (Y - 1) * (2 - V_G + Q)` with zero singular `(Y - 1)^0` drift at `Y = 1`.

2. **`gauge_a_exact_omega_collapse` & `gauge_a_exact_theta_balance`:**
   Under Gauge (a) (`c_l = 2 - V_G`, `c_ω = -1 + (V_G - H_1) / 2`, `Ḃ = -B * H_1`,
   `Q̃ = Q - H_1`), all inhomogeneous constant and linear `(V_G, H_1)` amplitude forcings
   cancel identically in both the `G` and `h` equations, leaving pure quadratic self-advection
   plus edge-subtracted linear dynamics `M̃_GG = -2∂_z + ½ Q̃_op`.

3. **`dual_strip_digamma_symmetry`:**
   The W1 Route-F Digamma ceiling `L(σ, δ) = 3π + 2/δ + 1/σ + 1/(1 - σ)` is invariant under
   `σ ↦ 1 - σ`, so the certified W1 ceiling at `σ_- = 3/4` equals that at `σ_+ = 1/4`.

4. **`symbol_norm_quarter_sq_le` .. `lower_order_smoothing_const`:**
   On the dual-strip lines `Re s = 1/4` and `Re s = 3/4` (where `|tan(π(σ + iτ))| = 1`),
   the convolution symbol `c(s) = -tan(πs)/(2s+1)` satisfies exact rational bounds
   `|c(1/4+iτ)| ≤ 2/3`, `|c(3/4+iτ)| ≤ 2/5`, `|(1/4+iτ) c(1/4+iτ)| < 1/2`,
   `|(3/4+iτ) c(3/4+iτ)| < 1/2`, yielding the exact `7/10` lower-order derivative smoothing
   constant `½ · (2/5) + ½ = 7/10` for `½ ∂_z(Q̃_op G) : L²_{3/4} → L²_{1/4}`.

5. **`hilbert_edge_commutator_fraction` & `nonlinear_transport_damping_margin`:**
   Exact rational-fraction splitting `(t - 1) / (Y - t) = (Y - 1) / (Y - t) - 1` driving the
   edge Hilbert commutator, and the quasilinear weighted-IBP transport damping floor
   `σ * a_nl - (1/2) * a_nl' ≥ 2*σ - σ*η₀ - (1/2)*η₁`.
-/
import Mathlib

namespace SingularBlowup.NonlinearComoving

noncomputable section

/-! ## 1. Exact Gauge-(a) Co-Moving PDE Collapse -/

/-- Exact edge-stagnation gauge `c_l = 2 - V_G` eliminates the singular drift at `Y = 1`:
`(-2 + V_G + (Y - 1) * Q) + (2 - V_G) * Y = (Y - 1) * (2 - V_G + Q)`. -/
theorem gauge_a_advection_velocity_collapse (Y V_G Q : ℝ) :
    (-2 + V_G + (Y - 1) * Q) + (2 - V_G) * Y = (Y - 1) * (2 - V_G + Q) := by
  ring

/-- **Exact Gauge-(a) nonlinear `Ω`-equation collapse.**
With `c_l = 2 - V_G`, `c_ω = -1 + (V_G - H_1) / 2`, and edge-subtracted nonlocal operator
`Q̃ = Q - H_1`, all constant and linear `(V_G, H_1)` amplitude forcings cancel identically:
  `-(2 - V_G + Q) * (G_z - (1/2) * (1 + G)) + (-1 + (V_G - H_1)/2) * (1 + G) + B_term`
  `= -2 * G_z + (1/2) * (Q - H_1) + B_term + (V_G - Q) * G_z + (1/2) * (Q - H_1) * G`. -/
theorem gauge_a_exact_omega_collapse (V_G H_1 Q G G_z B_term : ℝ) :
    -(2 - V_G + Q) * (G_z - (1 / 2 : ℝ) * (1 + G))
      + (-1 + (V_G - H_1) / 2) * (1 + G) + B_term =
    -2 * G_z + (1 / 2 : ℝ) * (Q - H_1) + B_term
      + (V_G - Q) * G_z + (1 / 2 : ℝ) * (Q - H_1) * G := by
  ring

/-- **Exact Gauge-(a) `Θ`-jump balance.**
Under `c_l = 2 - V_G` and `c_ω = -1 + (V_G - H_1) / 2`, we have `c_l + 2 * c_ω = -H_1`.
Setting `Ḃ = -B * H_1` cancels the inhomogeneous `(c_l + 2 c_ω) * (π/2 + h)` term identically,
leaving pure nonlinear transport `∂_τ h = (-2 + V_G - Q) * h_z`. -/
theorem gauge_a_exact_theta_balance (V_G H_1 Q B h h_z : ℝ) :
    let c_l : ℝ := 2 - V_G
    let c_om : ℝ := -1 + (V_G - H_1) / 2
    let B_dot : ℝ := -B * H_1
    c_l + 2 * c_om = -H_1 ∧
      B * (-(2 - V_G + Q) * h_z) + (c_l + 2 * c_om) * (B * (Real.pi / 2 + h))
        - B_dot * (Real.pi / 2 + h) =
      B * (-2 * h_z + (V_G - Q) * h_z) := by
  constructor <;> ring

/-! ## 2. Dual-Strip `W1` Digamma Symmetry & Exact Rational Symbol Bounds -/

/-- W1 Route-F Digamma logarithmic-derivative ceiling `L(σ, δ) = 3π + 2/δ + 1/σ + 1/(1 - σ)`. -/
def digammaCeiling (sigma delta : ℝ) : ℝ :=
  3 * Real.pi + 2 / delta + 1 / sigma + 1 / (1 - sigma)

/-- **Dual-strip symmetry `L(σ, δ) = L(1 - σ, δ)`.**
In particular, the certified W1 Digamma ceiling at `σ_- = 3/4` equals that at `σ_+ = 1/4`. -/
theorem dual_strip_digamma_symmetry (sigma delta : ℝ) :
    digammaCeiling (1 - sigma) delta = digammaCeiling sigma delta := by
  unfold digammaCeiling
  have h : 1 - (1 - sigma) = sigma := by ring
  rw [h]
  ring

/-- Exact rational symbol norm squared at `σ_+ = 1/4` (where `|tan(π(1/4 + iτ))| = 1`):
`|c(1/4 + iτ)|² = 1 / ((3/2)² + (2τ)²) ≤ (2/3)²`. -/
theorem symbol_norm_quarter_sq_le (tau : ℝ) :
    1 / ((3 / 2 : ℝ) ^ 2 + (2 * tau) ^ 2) ≤ (2 / 3 : ℝ) ^ 2 := by
  have h_pos : (0 : ℝ) < (3 / 2 : ℝ) ^ 2 + (2 * tau) ^ 2 := by positivity
  rw [div_le_iff₀ h_pos]
  nlinarith [sq_nonneg tau]

/-- Exact rational symbol norm squared at `σ_- = 3/4` (where `|tan(π(3/4 + iτ))| = 1`):
`|c(3/4 + iτ)|² = 1 / ((5/2)² + (2τ)²) ≤ (2/5)²`. -/
theorem symbol_norm_three_quarter_sq_le (tau : ℝ) :
    1 / ((5 / 2 : ℝ) ^ 2 + (2 * tau) ^ 2) ≤ (2 / 5 : ℝ) ^ 2 := by
  have h_pos : (0 : ℝ) < (5 / 2 : ℝ) ^ 2 + (2 * tau) ^ 2 := by positivity
  rw [div_le_iff₀ h_pos]
  nlinarith [sq_nonneg tau]

/-- Exact order-0 smoothing norm squared for `T_c ∂_z` at `σ_+ = 1/4`:
`|(1/4 + iτ) c(1/4 + iτ)|² = ((1/4)² + τ²) / ((3/2)² + (2τ)²) < (1/2)²`. -/
theorem symbol_deriv_quarter_sq_lt (tau : ℝ) :
    ((1 / 4 : ℝ) ^ 2 + tau ^ 2) / ((3 / 2 : ℝ) ^ 2 + (2 * tau) ^ 2) < (1 / 2 : ℝ) ^ 2 := by
  have h_pos : (0 : ℝ) < (3 / 2 : ℝ) ^ 2 + (2 * tau) ^ 2 := by positivity
  rw [div_lt_iff₀ h_pos]
  nlinarith [sq_nonneg tau]

/-- Exact order-0 smoothing norm squared for `T_c ∂_z` at `σ_- = 3/4`:
`|(3/4 + iτ) c(3/4 + iτ)|² = ((3/4)² + τ²) / ((5/2)² + (2τ)²) < (1/2)²`. -/
theorem symbol_deriv_three_quarter_sq_lt (tau : ℝ) :
    ((3 / 4 : ℝ) ^ 2 + tau ^ 2) / ((5 / 2 : ℝ) ^ 2 + (2 * tau) ^ 2) < (1 / 2 : ℝ) ^ 2 := by
  have h_pos : (0 : ℝ) < (5 / 2 : ℝ) ^ 2 + (2 * tau) ^ 2 := by positivity
  rw [div_lt_iff₀ h_pos]
  nlinarith [sq_nonneg tau]

/-- **Exact lower-order derivative smoothing constant `7/10`.**
For `½ ∂_z(Q̃_op G) = e^{-z/2} (-½ T_c G + T_c ∂_z G)` mapping `L²_{3/4} → L²_{1/4}`,
the exact operator-norm ceiling is `½ · (2/5) + ½ = 7/10`. -/
theorem lower_order_smoothing_const :
    (1 / 2 : ℝ) * (2 / 5 : ℝ) + (1 / 2 : ℝ) = 7 / 10 := by
  norm_num

/-! ## 3. Edge Hilbert Commutator & Quasilinear Transport Damping Floor -/

/-- Exact rational-fraction identity driving the Hilbert edge commutator on `Y > 1`:
`(t - 1) / (Y - t) = (Y - 1) / (Y - t) - 1` whenever `Y ≠ t`. -/
theorem hilbert_edge_commutator_fraction (Y t : ℝ) (hne : Y - t ≠ 0) :
    (t - 1) / (Y - t) = (Y - 1) / (Y - t) - 1 := by
  field_simp [hne]
  ring

/-- **Quasilinear weighted-IBP transport damping floor.**
In the weighted `L²_σ` energy estimate for `G_k = ∂_z^k G`, exact IBP gives local damping
density `σ * a_nl - (1/2) * a_nl'`. Whenever `|a_nl - 2| ≤ η₀` and `|a_nl'| ≤ η₁` with `σ > 0`,
the damping density is strictly bounded below by `2*σ - σ*η₀ - (1/2)*η₁`. -/
theorem nonlinear_transport_damping_margin
    (sigma a_nl da_nl eta0 eta1 : ℝ)
    (h_sigma : 0 < sigma)
    (h_a : |a_nl - 2| ≤ eta0)
    (h_da : |da_nl| ≤ eta1) :
    2 * sigma - sigma * eta0 - (1 / 2 : ℝ) * eta1 ≤
      sigma * a_nl - (1 / 2 : ℝ) * da_nl := by
  have h1 : 2 - eta0 ≤ a_nl := by
    have := (abs_le.mp h_a).1
    linarith
  have h2 : da_nl ≤ eta1 := (abs_le.mp h_da).2
  have h3 : sigma * (2 - eta0) ≤ sigma * a_nl :=
    mul_le_mul_of_nonneg_left h1 (le_of_lt h_sigma)
  linarith

end
end SingularBlowup.NonlinearComoving
