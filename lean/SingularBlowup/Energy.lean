/-
# Energy.lean — Track 1, lemmas 4–5 (weighted IBP + edge indicial roots)

The transport-dissipation identity (the exact integration-by-parts behind
`⟨G, -2∂_z G⟩_W = -∫ φ' W G² - [W G²]`, rigor_energy_estimate_setup.md §3) and
the edge indicial analysis (`findings_edge_condition_derivation.md` §5). Both are
self-contained 1D real analysis; no nonlocal / Hilbert machinery. They formalise
the *closed* skeleton of the Conjecture-2.4 proof; the nonlocal forcing coercivity
(★) is Track 2 and is deliberately NOT stated here.
-/
import Mathlib

namespace SingularBlowup

open Real MeasureTheory intervalIntegral

noncomputable section

/-! ## Lemma 4 — the weighted transport-dissipation identity (exact IBP) -/

/-- **Lemma 4.** For the transport operator `-2∂_z` in the weight `W = e^{-φ}`,
integration by parts is exact:
`  ∫_a^b ( e^{-φ}·G·(-2G') + φ'·e^{-φ}·G² ) dz = e^{-φ(a)}G(a)² - e^{-φ(b)}G(b)²,`
i.e. `⟨G, -2∂_z G⟩_W = -∫ φ' W G² - [W G²]_a^b`. The local dissipation rate `φ'`
appears on the left; the boundary term is the convective gap (`-2γ` for `φ=2γz`)
and the inflow piece vanishes under the edge/Robin condition `G(-∞)=0`. -/
theorem weighted_ibp (φ φ' G G' : ℝ → ℝ) (a b : ℝ)
    (hφ : ∀ z ∈ Set.uIcc a b, HasDerivAt φ (φ' z) z)
    (hG : ∀ z ∈ Set.uIcc a b, HasDerivAt G (G' z) z)
    (hcont : ContinuousOn
      (fun z => Real.exp (-φ z) * G z * (-2 * G' z)
        + φ' z * Real.exp (-φ z) * (G z) ^ 2) (Set.uIcc a b)) :
    (∫ z in a..b, Real.exp (-φ z) * G z * (-2 * G' z)
      + φ' z * Real.exp (-φ z) * (G z) ^ 2)
      = Real.exp (-φ a) * (G a) ^ 2 - Real.exp (-φ b) * (G b) ^ 2 := by
  have hderiv : ∀ z ∈ Set.uIcc a b,
      HasDerivAt (fun w => -(Real.exp (-φ w) * (G w) ^ 2))
        (Real.exp (-φ z) * G z * (-2 * G' z)
          + φ' z * Real.exp (-φ z) * (G z) ^ 2) z := by
    intro z hz
    have hd1 : HasDerivAt (fun w => Real.exp (-φ w)) (Real.exp (-φ z) * -(φ' z)) z :=
      ((hφ z hz).neg).exp
    have hd2 : HasDerivAt (fun w => (G w) ^ 2) (2 * (G z) ^ 1 * G' z) z :=
      (hG z hz).pow 2
    convert (hd1.mul hd2).neg using 1
    ring
  have hint : IntervalIntegrable
      (fun z => Real.exp (-φ z) * G z * (-2 * G' z)
        + φ' z * Real.exp (-φ z) * (G z) ^ 2) volume a b :=
    hcont.intervalIntegrable
  rw [integral_eq_sub_of_hasDerivAt hderiv hint]
  ring

/-! ## Lemma 5 — the edge indicial roots -/

/-- Edge power function `g_r(Y) = (Y-1)^r`. -/
def gpow (r Y : ℝ) : ℝ := (Y - 1) ^ r

theorem hasDerivAt_gpow (r : ℝ) {Y : ℝ} (hY : 1 < Y) :
    HasDerivAt (gpow r) (r * (Y - 1) ^ (r - 1)) Y := by
  have hY1 : (0 : ℝ) < Y - 1 := by linarith
  have h := ((hasDerivAt_id Y).sub_const 1).rpow_const (p := r) (Or.inl hY1.ne')
  simpa [gpow, one_mul, id_eq] using h

/-- **Lemma 5a (indicial reduction).** The homogeneous edge operator
`-2(Y-1)∂_Y - 1` (transport + `c_ω`) applied to `(Y-1)^r` returns the indicial
polynomial: `-2(Y-1)·g_r'(Y) - g_r(Y) = -(2r+1)·(Y-1)^r`. -/
theorem indicial_reduction (r : ℝ) {Y : ℝ} (hY : 1 < Y) :
    -2 * (Y - 1) * deriv (gpow r) Y - gpow r Y = -(2 * r + 1) * (Y - 1) ^ r := by
  have hY1 : (0 : ℝ) < Y - 1 := by linarith
  rw [(hasDerivAt_gpow r hY).deriv, gpow]
  have e : (Y - 1) * (Y - 1) ^ (r - 1) = (Y - 1) ^ r := by
    conv_rhs => rw [show r = 1 + (r - 1) by ring, Real.rpow_add hY1, Real.rpow_one]
  linear_combination (-2 * r) * e

/-- **Lemma 5b (indicial root).** `(Y-1)^r` is an edge eigenfunction at eigenvalue
`λ` (the indicial relation `-(2r+1) = λ`) iff `r = -(λ+1)/2`. (`lam` denotes `λ`;
for complex `λ` the same computation holds on the real exponent, with `r = Re`.) -/
theorem indicial_root_iff (r lam : ℝ) : -(2 * r + 1) = lam ↔ r = -(lam + 1) / 2 := by
  constructor <;> intro h <;> linarith

/-- **Lemma 5c (weighted admissibility).** The edge branch `(Y-1)^r` lies in the
weighted space (`G = (Y-1)^{1/2}·(Y-1)^r` decays, i.e. `r ≥ -1/2`) iff `λ ≤ 0`,
using `r = -(λ+1)/2`. This is the rigorous form of the Robin/edge condition that
drops the inadmissible homogeneous branch. -/
theorem indicial_admissible (lam : ℝ) :
    -(1 : ℝ) / 2 ≤ -(lam + 1) / 2 ↔ lam ≤ 0 := by
  constructor <;> intro h <;> linarith

end

end SingularBlowup
