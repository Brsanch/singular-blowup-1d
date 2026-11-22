/-
# ConjugationSymbol.lean — algebraic backbone of the W3 conjugation

Source: `docs/findings_session_2026_11_21_w3_conjugation.md`, ledger `LEM-10`.

The constrained linearization `M_σ = −2∂_z + K`, `(K g)^(s) = c(s+½) ĝ(s+½)`,
`c(s) = −tan(πs)/(2s+1)`, is conjugate to damped transport `A = −2∂_z − e^{−z/2}` through the
Mellin multiplier `N(s) = cos(πs) Γ(2s+2) / π`. This file checks the parts of that argument that
are identities about explicit functions:

1. `nSym_shift_mul_cSym` — the symbol identity `N(s+½) c(s+½) = −N(s)`.
2. `nSym_half` — `N(½) = 0` (the domain condition `ĝ(½) = 0`, and the reason `σ < ½`).
3. `flow_zero`, `flow_add` — the explicit damped-transport flow is a one-parameter semigroup.
4. `flow_along_characteristic`, `hasDerivAt_damp_along_characteristic` — along `z = z₀ + 2τ` the
   flow is `k(z₀)` times a factor whose τ-derivative is `−e^{−z/2}` times itself, i.e. the flow
   solves `∂_τ k = −2∂_z k − e^{−z/2} k` along characteristics.
5. `norm_flow_le` — the damping factor is at most 1 for `τ ≥ 0`.
6. `cRe_neg`, `cRe_pos` — on the real axis `c < 0` on `(0, ½)` and `c > 0` on `(½, 1)`, so
   `c(σ+¼) > 0` only for `σ > ¼` (correction to Lemma C of `convective_weight_spectrum.md`).
7. `reC_numerator_pos`, `reC_denominator_pos` — the sign factors of `Re c(x+it)` for `0 < x < ½`
   and every `t` (the closed formula for `Re c` itself is not formalized).

SCOPE. Nothing here is an operator statement. Not formalized: the two-sided bound
`|N(σ+it)| ≍ (1+|t|)^{2σ+3/2}`, the isomorphism `𝒩 : H^p_σ → L²_σ`, the domain identification
`D(M_σ) = 𝒩 D(A)`, the Sobolev multiplier bound for the damping factor, strong continuity, the
generator, and the constant `C*(σ)`. Those are on paper only (findings doc §§1–2).
-/
import Mathlib

namespace SingularBlowup.ConjugationSymbol

noncomputable section

/-- Ladder symbol `c(s) = −tan(πs)/(2s+1)`. -/
def cSym (s : ℂ) : ℂ := -Complex.tan ((Real.pi : ℂ) * s) / (2 * s + 1)

/-- Conjugating multiplier `N(s) = cos(πs) Γ(2s+2) / π`. -/
def nSym (s : ℂ) : ℂ := Complex.cos ((Real.pi : ℂ) * s) * Complex.Gamma (2 * s + 2) / (Real.pi : ℂ)

/-- **Symbol identity.** `N(s+½) c(s+½) = −N(s)` wherever `sin(πs) ≠ 0` and `2s+2 ≠ 0`. -/
theorem nSym_shift_mul_cSym (s : ℂ) (hsin : Complex.sin ((Real.pi : ℂ) * s) ≠ 0)
    (h2 : 2 * s + 2 ≠ 0) : nSym (s + 1 / 2) * cSym (s + 1 / 2) = -nSym s := by
  have hpi : (Real.pi : ℂ) ≠ 0 := Complex.ofReal_ne_zero.mpr Real.pi_ne_zero
  have harg : (Real.pi : ℂ) * (s + 1 / 2) = (Real.pi : ℂ) * s + (Real.pi : ℂ) / 2 := by ring
  have hc : Complex.cos ((Real.pi : ℂ) * (s + 1 / 2)) = -Complex.sin ((Real.pi : ℂ) * s) := by
    rw [harg, Complex.cos_add_pi_div_two]
  have hs : Complex.sin ((Real.pi : ℂ) * (s + 1 / 2)) = Complex.cos ((Real.pi : ℂ) * s) := by
    rw [harg, Complex.sin_add_pi_div_two]
  have hG : Complex.Gamma (2 * (s + 1 / 2) + 2) = (2 * s + 2) * Complex.Gamma (2 * s + 2) := by
    rw [show 2 * (s + 1 / 2) + 2 = (2 * s + 2) + 1 by ring, Complex.Gamma_add_one _ h2]
  have hden : 2 * (s + 1 / 2) + 1 = 2 * s + 2 := by ring
  have hs1 : s + 1 ≠ 0 := by
    intro h
    exact h2 (by linear_combination 2 * h)
  unfold nSym cSym
  rw [Complex.tan_eq_sin_div_cos, hc, hs, hG, hden]
  field_simp

/-- `N(½) = 0`: the zero that encodes the domain condition `ĝ(½) = 0`. -/
theorem nSym_half : nSym (1 / 2) = 0 := by
  have h : (Real.pi : ℂ) * (1 / 2) = (Real.pi : ℂ) / 2 := by ring
  unfold nSym
  rw [h, Complex.cos_pi_div_two]
  simp

/-- Damping factor `exp(−x e^{−z/2})`. -/
def damp (x z : ℝ) : ℝ := Real.exp (-(x * Real.exp (-z / 2)))

/-- Damped-transport flow `(T(τ)k)(z) = k(z−2τ) · exp(−(1−e^{−τ}) e^{−(z−2τ)/2})`. -/
def flow (τ : ℝ) (k : ℝ → ℂ) (z : ℝ) : ℂ :=
  k (z - 2 * τ) * ((damp (1 - Real.exp (-τ)) (z - 2 * τ) : ℝ) : ℂ)

theorem flow_zero (k : ℝ → ℂ) : flow 0 k = k := by
  funext z
  simp [flow, damp]

/-- Exponent bookkeeping behind the semigroup law. -/
theorem damp_mul_damp (s t w : ℝ) :
    damp (1 - Real.exp (-t)) w * damp (1 - Real.exp (-s)) (w + 2 * t)
      = damp (1 - Real.exp (-(s + t))) w := by
  have h1 : Real.exp (-(w + 2 * t) / 2) = Real.exp (-w / 2) * Real.exp (-t) := by
    rw [← Real.exp_add]; congr 1; ring
  have h2 : Real.exp (-(s + t)) = Real.exp (-s) * Real.exp (-t) := by
    rw [← Real.exp_add]; congr 1; ring
  unfold damp
  rw [← Real.exp_add, h1, h2]
  congr 1
  ring

/-- **Semigroup law** for the explicit flow. -/
theorem flow_add (s t : ℝ) (k : ℝ → ℂ) : flow s (flow t k) = flow (s + t) k := by
  funext z
  have hz : z - 2 * s - 2 * t = z - 2 * (s + t) := by ring
  have hw : z - 2 * s = (z - 2 * (s + t)) + 2 * t := by ring
  have key := damp_mul_damp s t (z - 2 * (s + t))
  simp only [flow]
  rw [hz, mul_assoc, ← Complex.ofReal_mul]
  rw [show damp (1 - Real.exp (-s)) (z - 2 * s)
        = damp (1 - Real.exp (-s)) ((z - 2 * (s + t)) + 2 * t) by rw [← hw]]
  rw [key]

/-- Along the characteristic `z = z₀ + 2τ` the flow is `k(z₀)` times the damping factor. -/
theorem flow_along_characteristic (τ z0 : ℝ) (k : ℝ → ℂ) :
    flow τ k (z0 + 2 * τ) = k z0 * ((damp (1 - Real.exp (-τ)) z0 : ℝ) : ℂ) := by
  simp [flow]

/-- The damping factor along a characteristic satisfies `d/dτ D = −e^{−(z₀+2τ)/2} · D`:
the flow solves `∂_τ k = −2∂_z k − e^{−z/2} k` along characteristics. -/
theorem hasDerivAt_damp_along_characteristic (τ z0 : ℝ) :
    HasDerivAt (fun τ' : ℝ => damp (1 - Real.exp (-τ')) z0)
      (-Real.exp (-(z0 + 2 * τ) / 2) * damp (1 - Real.exp (-τ)) z0) τ := by
  have hsplit : Real.exp (-(z0 + 2 * τ) / 2) = Real.exp (-z0 / 2) * Real.exp (-τ) := by
    rw [← Real.exp_add]; congr 1; ring
  have h1 : HasDerivAt (fun τ' : ℝ => Real.exp (-τ')) (Real.exp (-τ) * (-1)) τ :=
    (hasDerivAt_neg τ).exp
  have h2 : HasDerivAt (fun τ' : ℝ => -((1 - Real.exp (-τ')) * Real.exp (-z0 / 2)))
      (-((-(Real.exp (-τ) * (-1))) * Real.exp (-z0 / 2))) τ :=
    ((h1.const_sub 1).mul_const _).neg
  have h3 := h2.exp
  unfold damp
  convert h3 using 1
  rw [hsplit]
  ring

theorem damp_pos (x z : ℝ) : 0 < damp x z := Real.exp_pos _

theorem damp_le_one {x : ℝ} (hx : 0 ≤ x) (z : ℝ) : damp x z ≤ 1 := by
  unfold damp
  rw [Real.exp_le_one_iff]
  have := mul_nonneg hx (Real.exp_pos (-z / 2)).le
  linarith

/-- For `τ ≥ 0` the flow is pointwise dominated by the translated datum. -/
theorem norm_flow_le {τ : ℝ} (hτ : 0 ≤ τ) (k : ℝ → ℂ) (z : ℝ) :
    ‖flow τ k z‖ ≤ ‖k (z - 2 * τ)‖ := by
  have hx : 0 ≤ 1 - Real.exp (-τ) := by
    have : Real.exp (-τ) ≤ 1 := Real.exp_le_one_iff.mpr (by linarith)
    linarith
  have hd := damp_le_one hx (z - 2 * τ)
  have hp := damp_pos (1 - Real.exp (-τ)) (z - 2 * τ)
  unfold flow
  rw [norm_mul, Complex.norm_real, Real.norm_of_nonneg hp.le]
  exact mul_le_of_le_one_right (norm_nonneg _) hd

/-- The ladder symbol on the real axis. -/
def cRe (x : ℝ) : ℝ := -Real.tan (Real.pi * x) / (2 * x + 1)

theorem cSym_ofReal (x : ℝ) : cSym (x : ℂ) = ((cRe x : ℝ) : ℂ) := by
  simp [cSym, cRe]

/-- `c < 0` on `(0, ½)`. -/
theorem cRe_neg {x : ℝ} (h0 : 0 < x) (h1 : x < 1 / 2) : cRe x < 0 := by
  have hpi := Real.pi_pos
  have ht : 0 < Real.tan (Real.pi * x) :=
    Real.tan_pos_of_pos_of_lt_pi_div_two (by positivity) (by nlinarith)
  unfold cRe
  exact div_neg_of_neg_of_pos (by linarith) (by linarith)

/-- `c > 0` on `(½, 1)`. With `x = σ + ¼` this is `σ > ¼`: the positive band of Lemma C exists
only for `σ > ¼`. -/
theorem cRe_pos {x : ℝ} (h0 : 1 / 2 < x) (h1 : x < 1) : 0 < cRe x := by
  have hpi := Real.pi_pos
  have hper : Real.tan (Real.pi * x) = Real.tan (Real.pi * x - Real.pi) := by
    have := Real.tan_periodic (Real.pi * x - Real.pi)
    rw [sub_add_cancel] at this
    exact this
  have ht : Real.tan (Real.pi * x - Real.pi) < 0 :=
    Real.tan_neg_of_neg_of_pi_div_two_lt (by nlinarith) (by nlinarith)
  unfold cRe
  rw [hper]
  exact div_pos (by linarith) (by linarith)

/-! ### Sign of `Re c(x+it)` off the real axis

`Re c(x+it) = −[(2x+1) sin(2πx) + 2t sinh(2πt)] / [(cos 2πx + cosh 2πt) |2(x+it)+1|²]`
(the identity itself is checked numerically, probe WC-4, not here). The two lemmas below give
the sign of the numerator and of the first denominator factor for `0 < x < ½`, hence
`Re c(x+it) < 0` for every `t` there. With `x = σ + ¼` this is the dissipative range `σ < ¼`. -/

theorem mul_sinh_nonneg (t : ℝ) {a : ℝ} (ha : 0 < a) : 0 ≤ t * Real.sinh (a * t) := by
  rcases le_total 0 t with ht | ht
  · exact mul_nonneg ht (Real.sinh_nonneg_iff.mpr (mul_nonneg ha.le ht))
  · exact mul_nonneg_of_nonpos_of_nonpos ht
      (Real.sinh_nonpos_iff.mpr (mul_nonpos_of_nonneg_of_nonpos ha.le ht))

/-- Numerator of `−Re c(x+it)` is positive for `0 < x < ½`, every `t`. -/
theorem reC_numerator_pos {x : ℝ} (h0 : 0 < x) (h1 : x < 1 / 2) (t : ℝ) :
    0 < (2 * x + 1) * Real.sin (2 * Real.pi * x) + 2 * t * Real.sinh (2 * Real.pi * t) := by
  have hpi := Real.pi_pos
  have hsin : 0 < Real.sin (2 * Real.pi * x) :=
    Real.sin_pos_of_pos_of_lt_pi (by positivity) (by nlinarith)
  have hsh : 0 ≤ t * Real.sinh (2 * Real.pi * t) := mul_sinh_nonneg t (by positivity)
  have h1' : 0 < (2 * x + 1) * Real.sin (2 * Real.pi * x) := mul_pos (by linarith) hsin
  linarith

/-- `cos(2πx) + cosh(2πt) > 0` for `0 < x < ½`, every `t`. -/
theorem reC_denominator_pos {x : ℝ} (h0 : 0 < x) (h1 : x < 1 / 2) (t : ℝ) :
    0 < Real.cos (2 * Real.pi * x) + Real.cosh (2 * Real.pi * t) := by
  have hpi := Real.pi_pos
  have hcos : 0 < Real.cos (Real.pi * x) :=
    Real.cos_pos_of_mem_Ioo ⟨by nlinarith, by nlinarith⟩
  have h2 : Real.cos (2 * Real.pi * x) = 2 * Real.cos (Real.pi * x) ^ 2 - 1 := by
    rw [show 2 * Real.pi * x = 2 * (Real.pi * x) by ring, Real.cos_two_mul]
  have hch : 1 ≤ Real.cosh (2 * Real.pi * t) := Real.one_le_cosh _
  have hsq : 0 < Real.cos (Real.pi * x) ^ 2 := by positivity
  rw [h2]
  linarith

end

end SingularBlowup.ConjugationSymbol
