/-
# Profile.lean — Track 1, lemmas 1–3 (the steady profile + edge algebra)

Self-contained real-analysis facts about the singular Hou–Luo self-similar
profile (Thm 2.3 / eq. 2.4 of arXiv:2604.01868). No Hilbert / Sobolev machinery:
these are the *closed* algebraic parts of the Conjecture-2.4 stability proof.

Numerical counterparts (this repo):
  • Lemma 1  ↔ validate/verify_singular_steady.py  (steady residual ~2e-16)
  • Lemma 2  ↔ the U(1) = -2 / (c_l+2c_ω)=0 edge balance
  • Lemma 3  ↔ src/stability_edge_z.py  (transport+c_ω collapse, ~1e-11)
-/
import Mathlib

namespace SingularBlowup

open Real

noncomputable section

/-! ## Profile constants and the steady drift -/

/-- Self-similar scaling exponent `c_l = 2` (Thm 2.3). -/
def cl : ℝ := 2

/-- Self-similar scaling exponent `c_ω = -1` (Thm 2.3). -/
def cω : ℝ := -1

/-- Far-field value of the self-similar velocity: `Ū ≡ -2` on `X > 1`. -/
def UbarFar : ℝ := -2

/-- The singular vorticity profile on the smooth region `X > 1`:
`Ω̄(Y) = (Y-1)^(-1/2)`. -/
def Ωbar (Y : ℝ) : ℝ := (Y - 1) ^ (-(1 : ℝ) / 2)

/-- The density/temperature profile, constant `π/2` on `X > 1`. -/
def Θbar : ℝ := Real.pi / 2

/-- Steady drift coefficient on `X > 1`: `ā(Y) = Ū + c_l·Y = -2 + 2Y`. -/
def drift (Y : ℝ) : ℝ := UbarFar + cl * Y

theorem drift_eq (Y : ℝ) : drift Y = 2 * (Y - 1) := by
  simp only [drift, UbarFar, cl]; ring

/-! ## Lemma 1 — the steady residual vanishes on `Y > 1` -/

/-- `Ω̄'(Y) = -(1/2)(Y-1)^(-3/2)` on `Y > 1`. -/
theorem hasDerivAt_Ωbar {Y : ℝ} (hY : 1 < Y) :
    HasDerivAt Ωbar (-(1 : ℝ) / 2 * (Y - 1) ^ (-(3 : ℝ) / 2)) Y := by
  have hY1 : (0 : ℝ) < Y - 1 := by linarith
  have h := ((hasDerivAt_id Y).sub_const 1).rpow_const (p := -(1 : ℝ) / 2)
    (Or.inl hY1.ne')
  have e : -(1 : ℝ) / 2 - 1 = -(3 : ℝ) / 2 := by norm_num
  simpa [Ωbar, one_mul, e, id_eq] using h

/-- **Lemma 1 (Ω-equation).** The steady residual
`R_Ω = ā·Ω̄' - c_ω·Ω̄ - Θ̄'` of eq. (2.4) vanishes identically on `Y > 1`
(`Θ̄' = 0` there). This is the exact form of the `~2e-16` numerical gate. -/
theorem steady_residual_Ω {Y : ℝ} (hY : 1 < Y) :
    drift Y * (-(1 : ℝ) / 2 * (Y - 1) ^ (-(3 : ℝ) / 2)) - cω * Ωbar Y - 0 = 0 := by
  have hY1 : (0 : ℝ) < Y - 1 := by linarith
  have e1 : (Y - 1) * (Y - 1) ^ (-(3 : ℝ) / 2) = (Y - 1) ^ (-(1 : ℝ) / 2) := by
    rw [show (-(1 : ℝ) / 2) = 1 + (-(3 : ℝ) / 2) by norm_num, Real.rpow_add hY1,
      Real.rpow_one]
  rw [drift_eq, Ωbar, cω]
  have collapse : 2 * (Y - 1) * (-(1 : ℝ) / 2 * (Y - 1) ^ (-(3 : ℝ) / 2))
      = -((Y - 1) * (Y - 1) ^ (-(3 : ℝ) / 2)) := by ring
  rw [collapse, e1]; ring

/-- **Lemma 1 (Θ-equation).** The steady residual
`R_Θ = ā·Θ̄' - (c_l+2c_ω)·Θ̄` vanishes on `Y > 1`: both `Θ̄' = 0` and
`c_l + 2c_ω = 0`. -/
theorem steady_residual_Θ (Y : ℝ) :
    drift Y * 0 - (cl + 2 * cω) * Θbar = 0 := by
  simp only [cl, cω]; ring

/-! ## Lemma 2 — the edge / Θ-balance algebra -/

/-- The self-similar velocity `Ū`, piecewise: `2√(1-Y) - 2` on `Y < 1`, `-2`
otherwise. -/
def Ubar (Y : ℝ) : ℝ := if Y < 1 then 2 * Real.sqrt (1 - Y) - 2 else -2

/-- **Lemma 2a.** `Ū(1) = -2` (the far value; the left branch limit agrees,
see `Ubar_left_limit_at_one`). -/
theorem Ubar_one : Ubar 1 = -2 := by
  simp [Ubar]

/-- The `Y → 1⁻` branch `2√(1-Y) - 2` also evaluates to `-2` at `Y = 1`
(`√0 = 0`), so `Ū` is continuous across the edge. -/
theorem Ubar_left_limit_at_one : 2 * Real.sqrt (1 - 1) - 2 = -2 := by
  simp

/-- **Lemma 2b.** The Θ-balance constant vanishes: `c_l + 2c_ω = 0`.
This kills the RHS of the steady Θ-equation. -/
theorem theta_balance : cl + 2 * cω = 0 := by
  simp only [cl, cω]; ring

/-- **Lemma 2c.** The steady drift vanishes at the edge: `ā(1) = Ū(1) + c_l = 0`.
Together with `Θ̄' = (π/2)δ(X-1)` this makes the singular term
`ā·Θ̄'` have zero coefficient — the consistency condition that lets the singular
profile solve the distributional Θ-equation across the jump. -/
theorem drift_edge_zero : drift 1 = 0 := by
  simp only [drift, UbarFar, cl]; ring

/-- **Lemma 2 (consistency corollary).** The jump coefficient
`(Ū(1) + c_l)·Θ̄ = 0`. -/
theorem theta_jump_coeff_zero : (Ubar 1 + cl) * Θbar = 0 := by
  rw [Ubar_one]; simp only [cl]; ring

/-! ## Lemma 3 — the `c_ω`-absorption change of variables -/

/-- Physical shape `g(Y) = (Y-1)^(-1/2)·G(log(Y-1))` recovered from the
edge-resolved variable `G(z)`, `z = log(Y-1)`. -/
def gOf (G : ℝ → ℝ) (Y : ℝ) : ℝ := (Y - 1) ^ (-(1 : ℝ) / 2) * G (Real.log (Y - 1))

/-- Derivative of `g(Y) = (Y-1)^(-1/2)·G(log(Y-1))` (product + chain rule). -/
theorem hasDerivAt_gOf {G G' : ℝ → ℝ} (hG : ∀ z, HasDerivAt G (G' z) z)
    {Y : ℝ} (hY : 1 < Y) :
    HasDerivAt (gOf G)
      ((-(1 : ℝ) / 2 * (Y - 1) ^ (-(3 : ℝ) / 2)) * G (Real.log (Y - 1))
        + (Y - 1) ^ (-(1 : ℝ) / 2) * (G' (Real.log (Y - 1)) * (Y - 1)⁻¹)) Y := by
  have hY1 : (0 : ℝ) < Y - 1 := by linarith
  have hA : HasDerivAt (fun y => (y - 1) ^ (-(1 : ℝ) / 2))
      (-(1 : ℝ) / 2 * (Y - 1) ^ (-(3 : ℝ) / 2)) Y := by
    have h := ((hasDerivAt_id Y).sub_const 1).rpow_const (p := -(1 : ℝ) / 2)
      (Or.inl hY1.ne')
    have e : -(1 : ℝ) / 2 - 1 = -(3 : ℝ) / 2 := by norm_num
    simpa [one_mul, e, id_eq] using h
  have hlog : HasDerivAt (fun y => Real.log (y - 1)) ((Y - 1)⁻¹) Y := by
    have h := ((hasDerivAt_id Y).sub_const 1).log hY1.ne'
    simpa [one_div, id_eq] using h
  have hB : HasDerivAt (fun y => G (Real.log (y - 1)))
      (G' (Real.log (Y - 1)) * (Y - 1)⁻¹) Y :=
    (hG (Real.log (Y - 1))).comp Y hlog
  simpa [gOf] using hA.mul hB

/-- **Lemma 3 (`c_ω`-absorption).** In the edge-resolved variable the transport
`-ā g_Y` (with `ā = 2(Y-1)`) plus the source `c_ω g` (with `c_ω = -1`) collapses
to the clean `-2 ∂_z G`: the `+G` produced by `-ā g_Y` exactly cancels the `-G`
from `c_ω g`. Concretely, for `g(Y) = (Y-1)^(-1/2)·G(log(Y-1))`,
`  -2(Y-1)·g'(Y) - g(Y) = -2·(Y-1)^(-1/2)·G'(log(Y-1))`,
the RHS being `-2 ∂_z G` written back in `(Y-1)^(-1/2)` scale. -/
theorem cω_absorption {G G' : ℝ → ℝ} (hG : ∀ z, HasDerivAt G (G' z) z)
    {Y : ℝ} (hY : 1 < Y) :
    -2 * (Y - 1) * deriv (gOf G) Y - gOf G Y
      = -2 * (Y - 1) ^ (-(1 : ℝ) / 2) * G' (Real.log (Y - 1)) := by
  have hY1 : (0 : ℝ) < Y - 1 := by linarith
  rw [(hasDerivAt_gOf hG hY).deriv]
  simp only [gOf]
  have e1 : (Y - 1) * (Y - 1) ^ (-(3 : ℝ) / 2) = (Y - 1) ^ (-(1 : ℝ) / 2) := by
    rw [show (-(1 : ℝ) / 2) = 1 + (-(3 : ℝ) / 2) by norm_num, Real.rpow_add hY1,
      Real.rpow_one]
  have e2 : (Y - 1) * (Y - 1)⁻¹ = 1 := mul_inv_cancel₀ hY1.ne'
  linear_combination (G (Real.log (Y - 1))) * e1
    + (-2 * (Y - 1) ^ (-(1 : ℝ) / 2) * G' (Real.log (Y - 1))) * e2

end

end SingularBlowup
