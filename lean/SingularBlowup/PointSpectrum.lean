/-
# PointSpectrum.lean — LEM-08 integer obstruction (discrete half)

Formalizes the **exact integer incompatibility** at the core of LEM-08
(point-spectral stability mod 𝒩 of the frozen Ω-shape generator):

  line-L² of Ĝ_p · P forces k_max < 3/4 and |k_min| < 3/4
  ⟹ for integers k_max = k_min = 0 ⟹ Laurent degree d = 0;
  killing the s = 0 edge-constant pole requires d ≥ 1;
  d = 0 ∧ d ≥ 1 is impossible.

The continuous inputs (Stirling decay |Ĝ_p| ∼ e^{−(3π/2)|τ|}, residue Γ(λ) ≠ 0)
are machine-gated in Python (`validate/lem08_integer_obstruction.py`,
`validate/edge_pointspectrum_rigid.py`) and are **not** re-proved here — this
file is the firewall-clean discrete skeleton only.

Numerical / paper companions:
  • docs/ledger/LEM-08_point_spectrum_stability.md
  • docs/desk/edge_mellin_symbol.md §5.3–5.5
  • src/lem08_integer_obstruction.py
-/
import Mathlib

namespace SingularBlowup.PointSpectrum

open Int

noncomputable section

/-- The L²-line budget for Laurent modes: `k_max < 3/4` and `|k_min| < 3/4`,
coming from Stirling decay `3π/2` vs per-mode growth `2π`
(`3π/2 ÷ 2π = 3/4`). -/
def l2Budget : ℝ := (3 : ℝ) / 4

theorem l2Budget_eq : l2Budget = (3 : ℝ) / 4 := rfl

/-- **Integer L² admissibility.** Integers strictly inside `(-3/4, 3/4)` are
exactly `{0}`. -/
theorem int_lt_three_quarters_eq_zero (k : ℤ)
    (h : (k : ℝ) < (3 : ℝ) / 4)
    (h' : -((3 : ℝ) / 4) < (k : ℝ)) : k = 0 := by
  have hk1 : (k : ℝ) < (1 : ℝ) := by
    have : (3 : ℝ) / 4 < (1 : ℝ) := by norm_num
    linarith
  have hk2 : (-(1 : ℝ)) < (k : ℝ) := by
    have : (-(1 : ℝ)) < -((3 : ℝ) / 4) := by norm_num
    linarith
  have hlt : k < (1 : ℤ) := by
    exact_mod_cast hk1
  have hgt : (-1 : ℤ) < k := by
    exact_mod_cast hk2
  interval_cases k
  rfl

/-- **Only the constant Laurent mode is L²-admissible** among integers. -/
theorem only_constant_mode_l2_admissible (kMax kMin : ℤ)
    (hMax : (kMax : ℝ) < l2Budget)
    (hMin : |(kMin : ℝ)| < l2Budget)
    (hOrd : kMin ≤ kMax) :
    kMax = 0 ∧ kMin = 0 := by
  have hMax' : (kMax : ℝ) < (3 : ℝ) / 4 := by simpa [l2Budget] using hMax
  have hMinAbs : |(kMin : ℝ)| < (3 : ℝ) / 4 := by simpa [l2Budget] using hMin
  have hMinLo : -((3 : ℝ) / 4) < (kMin : ℝ) := neg_lt_of_abs_lt hMinAbs
  have hMinHi : (kMin : ℝ) < (3 : ℝ) / 4 := lt_of_abs_lt hMinAbs
  have hMaxLo : -((3 : ℝ) / 4) < (kMax : ℝ) := by
    have hkm : (kMin : ℝ) ≤ (kMax : ℝ) := Int.cast_le.mpr hOrd
    linarith
  have hkMax0 : kMax = 0 := int_lt_three_quarters_eq_zero kMax hMax' hMaxLo
  have hkMin0 : kMin = 0 := int_lt_three_quarters_eq_zero kMin hMinHi hMinLo
  exact ⟨hkMax0, hkMin0⟩

/-- Laurent degree of a two-sided integer mode support. -/
def laurentDegree (kMax kMin : ℤ) : ℤ := kMax - kMin

theorem laurentDegree_nonneg (kMax kMin : ℤ) (h : kMin ≤ kMax) :
    0 ≤ laurentDegree kMax kMin := by
  simpa [laurentDegree] using sub_nonneg.mpr h

/-- **Discrete LEM-08 obstruction.** Line-L² forces the constant mode
(`d = 0`); edge-const kill forces `d ≥ 1`; incompatible. -/
theorem integer_obstruction
    (kMax kMin : ℤ)
    (hMax : (kMax : ℝ) < l2Budget)
    (hMin : |(kMin : ℝ)| < l2Budget)
    (hOrd : kMin ≤ kMax)
    (hEdgeKill : (1 : ℤ) ≤ laurentDegree kMax kMin) :
    False := by
  rcases only_constant_mode_l2_admissible kMax kMin hMax hMin hOrd with ⟨hM, hN⟩
  subst hM; subst hN
  simp [laurentDegree] at hEdgeKill

/-- Specialization: the constant mode cannot also kill the edge pole. -/
theorem constant_mode_cannot_kill_edge :
    ¬ ((1 : ℤ) ≤ laurentDegree (0 : ℤ) 0) := by
  intro h
  simp [laurentDegree] at h

/-- **Corollary (point-spectrum discrete half).** There is no integer Laurent
support that is simultaneously L²-admissible and edge-pole-killing. -/
theorem no_l2_and_edge_kill (kMax kMin : ℤ)
    (hMax : (kMax : ℝ) < l2Budget)
    (hMin : |(kMin : ℝ)| < l2Budget)
    (hOrd : kMin ≤ kMax) :
    ¬ ((1 : ℤ) ≤ laurentDegree kMax kMin) := by
  intro hKill
  exact integer_obstruction kMax kMin hMax hMin hOrd hKill

end
end SingularBlowup.PointSpectrum
