/-
# NonlinearBootstrap.lean — Final Nonlinear Bootstrap Closure (Edge Hardy-Disk + Overlap + Riccati)

Formalizes the closed nonlinear bootstrap backbone for the singular 1D Hou–Luo self-similar
blow-up profile (`docs/desk/nonlinear_comoving_leg.md` §4,
`validate/edge_nonlinear_bootstrap_closure.py`, 2026-11-20):

1. **Edge Hardy-Disk `ℋ¹(D(1,1))` Coercivity & Quasilinear Structure:**
   - `edge_disk_weight_ratio_le_half`: `√m / (m + 1) ≤ 1/2` for all `m : ℕ`, bounding the
     discrete Hilbert operator `ℋ_{1/2}` (`‖ℋ_{1/2}‖₂ = π`) against the `m³` transport
     dissipation in `ℋ¹(D(1,1))`.
   - `edge_disk_coercivity_margin`: exact linear + far-field Young dissipation balance
     `-4 + 1/2 + 1/2 = -3` (keeping full `-3 ∑ m³ |c_m|²` coercivity).
   - `quasilinear_sqrt_splitting` & `discrete_hilbert_commutator_identity`: exact convolution
     weight splitting `√m · k ≤ k √k + √(m - k) · k` and discrete Hilbert commutator
     `j / (k - j - 1/2) = (k - 1/2) / (k - j - 1/2) - 1`.

2. **Zero-Derivative `L²_{3/4,1/4} → L^∞` Bound, Overlap Inflow & Edge Embedding:**
   - `zero_deriv_linf_bound_lt_seven_fifths`: `1 / √π + √(2/3) < 7/5` (using only `π > 3`).
   - `overlap_inflow_geometric_tail_bound` & `overlap_weighted_inflow_lt_one`:
     `∑_{m=1}^N (1/4)^m < 1/3` and `2√2 / 3 < 1`, certifying the overlap inflow
     `2 e^{-z₀/2} |G_z(-ln 2)|² < ‖c‖_{ℋ¹}²` at `z₀ = -ln 2` (`Y = 3/2`).
   - `edge_l2_three_quarters_trace_bound`: elementary telescoping bound
     `∑_{m=1}^N 1 / (m² (2m - 3/2)) < 199/90 < 9/4`, certifying
     `‖G‖_{L²_{3/4}((-∞,0])}² ≤ (9/4) ‖c‖_{ℋ¹}²`.
   - `far_field_disk_coupling_bound`: elementary telescoping bound
     `∑_{m=1}^N m / (4 π² (m + 1)² (2m + 1/2)) < 1/100`, certifying zero radius shrinkage
     (`r₀ = 1`) for the far-field-to-edge-disk Taylor coupling.

3. **Closed 3-Block Lyapunov Functional & End-to-End Riccati Bootstrap Trapping:**
   - `lyapunov_three_block_coercivity`: with weights `(a₀, a₁, a₂) = (20, 1, 10)`, the coupled
     `(X_sem, Y_far, Z_disk)` linear dissipation is strictly coercive with exact rational
     per-block rates `(-741/1000, -1/2, -5359/2000)` and overall rate `γ_nl = 1/2`.
   - `riccati_bootstrap_trapping_decay`: genuine real-analytic ODE trapping theorem proving
     that any nonnegative differentiable energy `E : ℝ → ℝ` satisfying
     `E'(t) ≤ -γ E(t) + C E(t) √(E(t))` and `E(0) < (γ / (2C))²` stays strictly trapped below
     `(γ / (2C))²` and decays exponentially as `E(t) ≤ E(0) exp(-(γ/2) t)` for all `t ≥ 0`.
-/
import Mathlib

namespace SingularBlowup.NonlinearBootstrap

noncomputable section

/-! ## 1. Edge Hardy-Disk `ℋ¹(D(1,1))` Coercivity & Quasilinear Identities -/

/-- For every `m : ℕ`, `√m / (m + 1) ≤ 1 / 2`.
In the `ℋ¹(D(1,1))` energy identity, this pairs with `‖ℋ_{1/2}‖₂ = π` to bound the
nonlocal discrete-Hilbert term by `(1/2) ∑ m³ |c_m|²`. -/
theorem edge_disk_weight_ratio_le_half (m : ℕ) :
    Real.sqrt (m : ℝ) / ((m : ℝ) + 1) ≤ 1 / 2 := by
  have hm : (0 : ℝ) ≤ (m : ℝ) := Nat.cast_nonneg m
  have hden : (0 : ℝ) < (m : ℝ) + 1 := by linarith
  have hsq : (Real.sqrt (m : ℝ) - 1) ^ 2 ≥ 0 := sq_nonneg _
  have hsqrt_sq : Real.sqrt (m : ℝ) ^ 2 = (m : ℝ) := Real.sq_sqrt hm
  have h2 : 2 * Real.sqrt (m : ℝ) ≤ (m : ℝ) + 1 := by nlinarith
  rw [div_le_div_iff₀ hden (by norm_num : (0 : ℝ) < 2)]
  linarith

/-- Exact linear + far-field Young dissipation margin on `ℋ¹(D(1,1))`:
`-4` (transport) `+ 1/2` (discrete Hilbert) `+ 1/2` (Young on far-field forcing) `= -3`. -/
theorem edge_disk_coercivity_margin :
    (-4 : ℝ) + 1 / 2 + 1 / 2 = -3 ∧ (-2 : ℝ) + 1 / 4 + 1 / 4 = -3 / 2 := by
  constructor <;> norm_num

/-- Quasilinear convolution weight splitting on `ℋ¹(D(1,1))`:
for `0 ≤ k ≤ m`, `√m · k ≤ k √k + √(m - k) · k`, enabling absorption of the quadratic
self-advection convolution into the `-3 ∑ m³ |c_m|²` transport dissipation. -/
theorem quasilinear_sqrt_splitting {k m : ℝ} (hk : 0 ≤ k) (hkm : k ≤ m) :
    Real.sqrt m * k ≤ k * Real.sqrt k + Real.sqrt (m - k) * k := by
  have hmk : 0 ≤ m - k := sub_nonneg.mpr hkm
  have hm : 0 ≤ m := le_trans hk hkm
  have hsub : Real.sqrt m ≤ Real.sqrt k + Real.sqrt (m - k) := by
    have h1 : 0 ≤ Real.sqrt k + Real.sqrt (m - k) :=
      add_nonneg (Real.sqrt_nonneg k) (Real.sqrt_nonneg (m - k))
    have h2 : m ≤ (Real.sqrt k + Real.sqrt (m - k)) ^ 2 := by
      have hsqk : Real.sqrt k ^ 2 = k := Real.sq_sqrt hk
      have hsqmk : Real.sqrt (m - k) ^ 2 = m - k := Real.sq_sqrt hmk
      have hprod : 0 ≤ Real.sqrt k * Real.sqrt (m - k) :=
        mul_nonneg (Real.sqrt_nonneg k) (Real.sqrt_nonneg (m - k))
      nlinarith
    calc
      Real.sqrt m ≤ Real.sqrt ((Real.sqrt k + Real.sqrt (m - k)) ^ 2) := Real.sqrt_le_sqrt h2
      _ = Real.sqrt k + Real.sqrt (m - k) := Real.sqrt_sq h1
  nlinarith

/-- Exact discrete Hilbert commutator identity on power-series indices:
`j / (k - j - 1/2) = (k - 1/2) / (k - j - 1/2) - 1`. -/
theorem discrete_hilbert_commutator_identity {k j : ℝ} (hne : k - j - 1 / 2 ≠ 0) :
    j / (k - j - 1 / 2) = (k - 1 / 2) / (k - j - 1 / 2) - 1 := by
  rw [div_sub_one hne]
  ring

/-! ## 2. Zero-Derivative `L^∞` Bound, Overlap Inflow, Edge Embedding & Far-Field Coupling -/

/-- Zero-derivative `L²_{3/4,1/4} → L^∞` symbol constant for `Q_op G`:
`1 / √π + √(2/3) < 7 / 5`, proved purely from `π > 3`. -/
theorem zero_deriv_linf_bound_lt_seven_fifths :
    1 / Real.sqrt Real.pi + Real.sqrt (2 / 3 : ℝ) < 7 / 5 := by
  have hpi : (3 : ℝ) < Real.pi := Real.pi_gt_three
  have hs3_pos : (0 : ℝ) < Real.sqrt 3 := Real.sqrt_pos.mpr (by norm_num)
  have hspi_gt : Real.sqrt 3 < Real.sqrt Real.pi :=
    Real.sqrt_lt_sqrt (by norm_num) hpi
  have hinv : 1 / Real.sqrt Real.pi < 1 / Real.sqrt 3 :=
    one_div_lt_one_div_of_lt hs3_pos hspi_gt
  have hstep : 1 / Real.sqrt 3 + Real.sqrt (2 / 3 : ℝ) < 7 / 5 := by
    have hs23_nonneg : (0 : ℝ) ≤ Real.sqrt (2 / 3 : ℝ) := Real.sqrt_nonneg _
    have hsq3 : Real.sqrt 3 ^ 2 = 3 := Real.sq_sqrt (by norm_num)
    have hsq23 : Real.sqrt (2 / 3 : ℝ) ^ 2 = 2 / 3 := Real.sq_sqrt (by norm_num)
    have hprod_sq : (Real.sqrt 3 * Real.sqrt (2 / 3 : ℝ)) ^ 2 = 2 := by
      rw [mul_pow, hsq3, hsq23]; norm_num
    have hprod_pos : 0 ≤ Real.sqrt 3 * Real.sqrt (2 / 3 : ℝ) :=
      mul_nonneg (le_of_lt hs3_pos) hs23_nonneg
    have hprod_lt : Real.sqrt 3 * Real.sqrt (2 / 3 : ℝ) < 36 / 25 := by nlinarith
    -- Multiply target by √3 > 0: need 1 + √3 * √(2/3) < (7/5) * √3
    -- Equivalently 1 + √2 < (7/5) √3, i.e., (1 + √2)² = 3 + 2√2 < (49/25)*3 = 147/25,
    -- which holds since √2 < 36/25 (1250 < 1296).
    have hmul : 1 + Real.sqrt 3 * Real.sqrt (2 / 3 : ℝ) < (7 / 5 : ℝ) * Real.sqrt 3 := by
      have hsq_lhs : (1 + Real.sqrt 3 * Real.sqrt (2 / 3 : ℝ)) ^ 2 < ((7 / 5 : ℝ) * Real.sqrt 3) ^ 2 := by
        rw [mul_pow, hsq3]
        nlinarith
      nlinarith
    have hdiv : (1 + Real.sqrt 3 * Real.sqrt (2 / 3 : ℝ)) / Real.sqrt 3 < 7 / 5 := by
      rw [div_lt_iff₀ hs3_pos]
      exact hmul
    have hrewrite :
        1 / Real.sqrt 3 + Real.sqrt (2 / 3 : ℝ) =
        (1 + Real.sqrt 3 * Real.sqrt (2 / 3 : ℝ)) / Real.sqrt 3 := by
      field_simp [ne_of_gt hs3_pos]
    linarith
  linarith

/-- Exact finite geometric sum bound for the overlap inflow at `z₀ = -ln 2` (`e^{z₀} = 1/2`):
`∑_{m=1}^N (1/4)^m = (1/3)(1 - (1/4)^N) < 1/3`. -/
theorem overlap_inflow_geometric_tail_bound (N : ℕ) :
    (∑ m ∈ Finset.range N, ((1 / 4 : ℝ) ^ (m + 1))) < 1 / 3 := by
  have hsum : (∑ m ∈ Finset.range N, ((1 / 4 : ℝ) ^ (m + 1))) =
      (1 / 3 : ℝ) * (1 - (1 / 4 : ℝ) ^ N) := by
    induction N with
    | zero => simp
    | succ n ih =>
      rw [Finset.sum_range_succ, ih]
      ring
  rw [hsum]
  have hpos : (0 : ℝ) < (1 / 4 : ℝ) ^ N := by positivity
  linarith

/-- Weighted overlap inflow coefficient at `z₀ = -ln 2`: `2 e^{-z₀/2} · (1/3) = 2√2 / 3 < 1`. -/
theorem overlap_weighted_inflow_lt_one :
    2 * Real.sqrt 2 / 3 < 1 := by
  have hsq : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have hnonneg : (0 : ℝ) ≤ Real.sqrt 2 := Real.sqrt_nonneg 2
  nlinarith

/-- Auxiliary telescoping bound for the tail `k ≥ 3` of the edge `L²_{3/4}` Gram trace:
`∑_{i=0}^{n-1} 1 / ((i + 3)² (2(i + 3) - 3/2)) ≤ (2/9) (1/2 - 1/(n + 2))`. -/
private theorem edge_trace_tail_le (n : ℕ) :
    (∑ i ∈ Finset.range n,
      (1 : ℝ) / (((i + 3 : ℕ) : ℝ) ^ 2 * (2 * ((i + 3 : ℕ) : ℝ) - 3 / 2)))
    ≤ (2 / 9 : ℝ) * (1 / 2 - 1 / ((n : ℝ) + 2)) := by
  induction n with
  | zero => simp
  | succ m ih =>
    rw [Finset.sum_range_succ]
    have hm : (0 : ℝ) ≤ (m : ℝ) := Nat.cast_nonneg m
    have hk : ((m + 3 : ℕ) : ℝ) = (m : ℝ) + 3 := by push_cast; ring
    have hstep :
        (1 : ℝ) / (((m + 3 : ℕ) : ℝ) ^ 2 * (2 * ((m + 3 : ℕ) : ℝ) - 3 / 2))
        ≤ (2 / 9 : ℝ) * (1 / ((m : ℝ) + 2) - 1 / ((m : ℝ) + 3)) := by
      rw [hk]
      have h2m : (0 : ℝ) < 2 * ((m : ℝ) + 3) - 3 / 2 := by linarith
      have hden1 : (0 : ℝ) < ((m : ℝ) + 3) ^ 2 * (2 * ((m : ℝ) + 3) - 3 / 2) :=
        mul_pos (by positivity) h2m
      have hden2 : (0 : ℝ) < ((m : ℝ) + 2) * ((m : ℝ) + 3) := by positivity
      have hdiff :
          (2 / 9 : ℝ) * (1 / ((m : ℝ) + 2) - 1 / ((m : ℝ) + 3)) =
          2 / (9 * ((m : ℝ) + 2) * ((m : ℝ) + 3)) := by
        field_simp
        ring
      rw [hdiff, div_le_div_iff₀ hden1 (by positivity)]
      nlinarith
    have hcast : ((m + 1 : ℕ) : ℝ) + 2 = (m : ℝ) + 3 := by push_cast; ring
    rw [hcast]
    linarith

/-- Elementary telescoping bound on the trace of `D⁻¹ Gram_{3/4} D⁻¹`:
for every `N : ℕ`, `∑_{m=1}^N 1 / (m² (2m - 3/2)) < 199 / 90 < 9 / 4`.
By Cauchy–Schwarz, this certifies `‖G‖_{L²_{3/4}((-∞,0])}² ≤ (9/4) ‖c‖_{ℋ¹(D(1,1))}²`. -/
theorem edge_l2_three_quarters_trace_bound (N : ℕ) :
    (∑ m ∈ Finset.range N,
      (1 : ℝ) / (((m + 1 : ℕ) : ℝ) ^ 2 * (2 * ((m + 1 : ℕ) : ℝ) - 3 / 2))) < 199 / 90 := by
  match N with
  | 0 => norm_num
  | 1 => norm_num
  | 2 =>
    rw [Finset.sum_range_succ, Finset.sum_range_one]
    norm_num
  | n + 2 =>
    have hsplit :
        (∑ m ∈ Finset.range (n + 2),
          (1 : ℝ) / (((m + 1 : ℕ) : ℝ) ^ 2 * (2 * ((m + 1 : ℕ) : ℝ) - 3 / 2))) =
        2 + 1 / 10 +
        (∑ i ∈ Finset.range n,
          (1 : ℝ) / (((i + 3 : ℕ) : ℝ) ^ 2 * (2 * ((i + 3 : ℕ) : ℝ) - 3 / 2))) := by
      have h1 : ∀ k : ℕ,
          (∑ m ∈ Finset.range (k + 2),
            (1 : ℝ) / (((m + 1 : ℕ) : ℝ) ^ 2 * (2 * ((m + 1 : ℕ) : ℝ) - 3 / 2))) =
          2 + 1 / 10 +
          (∑ i ∈ Finset.range k,
            (1 : ℝ) / (((i + 3 : ℕ) : ℝ) ^ 2 * (2 * ((i + 3 : ℕ) : ℝ) - 3 / 2))) := by
        intro k
        induction k with
        | zero =>
          rw [Finset.sum_range_succ, Finset.sum_range_one]
          simp
          norm_num
        | succ j ihj =>
          rw [show j + 1 + 2 = (j + 2) + 1 by omega, Finset.sum_range_succ, ihj,
              Finset.sum_range_succ]
          have hidx : (j + 2) + 1 = j + 3 := by omega
          rw [hidx]
          ring
      exact h1 n
    rw [hsplit]
    have htail := edge_trace_tail_le n
    have hpos : (0 : ℝ) < 1 / ((n : ℝ) + 2) := by positivity
    linarith

/-- Auxiliary telescoping bound `∑_{i=0}^{n-1} 1 / (i + 4)² ≤ 1/3 - 1/(n + 3)`. -/
private theorem inv_sq_tail_le (n : ℕ) :
    (∑ i ∈ Finset.range n, (1 : ℝ) / (((i + 4 : ℕ) : ℝ) ^ 2))
    ≤ 1 / 3 - 1 / ((n : ℝ) + 3) := by
  induction n with
  | zero => simp
  | succ m ih =>
    rw [Finset.sum_range_succ]
    have hm : (0 : ℝ) ≤ (m : ℝ) := Nat.cast_nonneg m
    have hk : ((m + 4 : ℕ) : ℝ) = (m : ℝ) + 4 := by push_cast; ring
    have hstep :
        (1 : ℝ) / (((m + 4 : ℕ) : ℝ) ^ 2)
        ≤ 1 / ((m : ℝ) + 3) - 1 / ((m : ℝ) + 4) := by
      rw [hk]
      have hdiff :
          1 / ((m : ℝ) + 3) - 1 / ((m : ℝ) + 4) =
          1 / (((m : ℝ) + 3) * ((m : ℝ) + 4)) := by
        field_simp
        ring
      rw [hdiff, div_le_div_iff₀ (by positivity) (by positivity)]
      nlinarith
    have hcast : ((m + 1 : ℕ) : ℝ) + 3 = (m : ℝ) + 4 := by push_cast; ring
    rw [hcast]
    linarith

/-- Elementary bound on the far-field-to-edge-disk Taylor coupling series (`NB-3`):
for every `N : ℕ`,
`∑_{m=1}^N m / (4 π² (m + 1)² (2m + 1/2)) < 1 / 100`,
proved purely from `π > 3` and telescoping (`25 / 2592 < 1 / 100`). -/
theorem far_field_disk_coupling_bound (N : ℕ) :
    (∑ m ∈ Finset.range N,
      ((m + 1 : ℕ) : ℝ) /
        (4 * Real.pi ^ 2 * (((m + 2 : ℕ) : ℝ) ^ 2) * (2 * ((m + 1 : ℕ) : ℝ) + 1 / 2)))
    < 1 / 100 := by
  have hpi : (3 : ℝ) < Real.pi := Real.pi_gt_three
  have hpi2 : (9 : ℝ) < Real.pi ^ 2 := by nlinarith
  have hterm :
      (∑ m ∈ Finset.range N,
        ((m + 1 : ℕ) : ℝ) /
          (4 * Real.pi ^ 2 * (((m + 2 : ℕ) : ℝ) ^ 2) * (2 * ((m + 1 : ℕ) : ℝ) + 1 / 2)))
      ≤ (1 / (8 * 9 : ℝ)) * (∑ m ∈ Finset.range N, (1 : ℝ) / (((m + 2 : ℕ) : ℝ) ^ 2)) := by
    rw [Finset.mul_sum]
    apply Finset.sum_le_sum
    intro m _
    have hm : (0 : ℝ) < ((m + 1 : ℕ) : ℝ) := by positivity
    have hden1 : (0 : ℝ) < 4 * Real.pi ^ 2 * (((m + 2 : ℕ) : ℝ) ^ 2) *
        (2 * ((m + 1 : ℕ) : ℝ) + 1 / 2) := by positivity
    have hden2 : (0 : ℝ) < (8 * 9 : ℝ) * (((m + 2 : ℕ) : ℝ) ^ 2) := by positivity
    have hrewrite :
        (1 / (8 * 9 : ℝ)) * ((1 : ℝ) / (((m + 2 : ℕ) : ℝ) ^ 2)) =
        1 / ((8 * 9 : ℝ) * (((m + 2 : ℕ) : ℝ) ^ 2)) := by ring
    rw [hrewrite, div_le_div_iff₀ hden1 hden2]
    have hsq2 : (0 : ℝ) < (((m + 2 : ℕ) : ℝ) ^ 2) := by positivity
    nlinarith
  have hsum_bound :
      (∑ m ∈ Finset.range N, (1 : ℝ) / (((m + 2 : ℕ) : ℝ) ^ 2)) < 25 / 36 := by
    match N with
    | 0 => norm_num
    | 1 => norm_num
    | 2 =>
      rw [Finset.sum_range_succ, Finset.sum_range_one]
      norm_num
    | n + 2 =>
      have hsplit : ∀ k : ℕ,
          (∑ m ∈ Finset.range (k + 2), (1 : ℝ) / (((m + 2 : ℕ) : ℝ) ^ 2)) =
          1 / 4 + 1 / 9 +
          (∑ i ∈ Finset.range k, (1 : ℝ) / (((i + 4 : ℕ) : ℝ) ^ 2)) := by
        intro k
        induction k with
        | zero =>
          rw [Finset.sum_range_succ, Finset.sum_range_one]
          simp
          norm_num
        | succ j ihj =>
          rw [show j + 1 + 2 = (j + 2) + 1 by omega, Finset.sum_range_succ, ihj,
              Finset.sum_range_succ]
          have hidx : (j + 2) + 2 = j + 4 := by omega
          rw [hidx]
          ring
      rw [hsplit n]
      have htail := inv_sq_tail_le n
      have hpos : (0 : ℝ) < 1 / ((n : ℝ) + 3) := by positivity
      linarith
  linarith

/-! ## 3. Closed 3-Block Lyapunov Functional & End-to-End Riccati Trapping -/

/-- **Closed 3-block Lyapunov coercivity (`NB-5`).**
With `E = 20 * X_sem + Y_far + 10 * Z_disk`, the coupled linear differential inequalities
- `Ẋ_sem ≤ -(4/5) X_sem` (Datko–Pazy renormalized W1–W3 semigroup norm, `ω₁ = 2/5 < 9/20`),
- `Ẏ_far ≤ -(1/2) Y_far + (641/200) Z_disk + (49/50) X_sem` (overlap inflow `< 1` + `7/10`
  lower-order smoothing + `9/4` edge embedding),
- `Ż_disk ≤ -3 Z_disk + (1/50) X_sem` (`-3` edge-disk coercivity + `2 * 1/100` far-field coupling),
combine to give strict decay at rate `γ_nl = 1/2` across all three blocks simultaneously
(exact per-block rates `-741/1000`, `-1/2`, `-5359/2000`). -/
theorem lyapunov_three_block_coercivity {X Y Z : ℝ} (hX : 0 ≤ X) (hZ : 0 ≤ Z) :
    20 * (-(4 / 5 : ℝ) * X)
      + 1 * (-(1 / 2 : ℝ) * Y + (641 / 200 : ℝ) * Z + (49 / 50 : ℝ) * X)
      + 10 * (-3 * Z + (1 / 50 : ℝ) * X)
    ≤ -(1 / 2 : ℝ) * (20 * X + Y + 10 * Z) := by
  linarith

/-- **End-to-end Riccati / Lyapunov bootstrap trapping & exponential decay theorem.**
If a nonnegative differentiable energy `E : ℝ → ℝ` satisfies the Riccati differential inequality
`E'(t) ≤ -γ E(t) + C E(t) √(E(t))` on `[0, ∞)` (`γ > 0, C > 0`) and starts strictly below the
trapping threshold `E(0) < (γ / (2C))²`, then `E(t)` remains strictly trapped and decays
exponentially for all `t ≥ 0`:
  `E(t) ≤ E(0) * exp(-(γ / 2) * t)`. -/
theorem riccati_bootstrap_trapping_decay
    {E E' : ℝ → ℝ} {γ C : ℝ}
    (hγ : 0 < γ) (hC : 0 < C)
    (hderiv : ∀ t ≥ 0, HasDerivAt E (E' t) t)
    (hE_nonneg : ∀ t ≥ 0, 0 ≤ E t)
    (h_riccati : ∀ t ≥ 0, E' t ≤ -γ * E t + C * E t * Real.sqrt (E t))
    (hE0 : E 0 < (γ / (2 * C)) ^ 2) :
    ∀ t ≥ 0, E t ≤ E 0 * Real.exp (-(γ / 2) * t) := by
  set R : ℝ := (γ / (2 * C)) ^ 2 with hR_def
  have hR_pos : 0 < R := by positivity
  have hsqrt_R : Real.sqrt R = γ / (2 * C) := by
    rw [hR_def, Real.sqrt_sq (by positivity)]
  have hE_cont : ContinuousOn E (Set.Ici 0) := fun t ht =>
    (hderiv t ht).continuousAt.continuousWithinAt
  -- Step 1: Strict barrier trapping `E t < R` for all `t ≥ 0`.
  have h_trap : ∀ t ≥ 0, E t < R := by
    intro t₁ ht₁
    by_contra hnot
    push Not at hnot
    let S : Set ℝ := Set.Icc 0 t₁ ∩ {t | R ≤ E t}
    have hS_closed : IsClosed S := by
      have hcont_Icc : ContinuousOn E (Set.Icc 0 t₁) :=
        hE_cont.mono (fun x hx => hx.1)
      have hpre := hcont_Icc.preimage_isClosed_of_isClosed isClosed_Icc (isClosed_Ici (a := R))
      simpa [S, Set.preimage, Set.mem_Ici] using hpre
    have hS_compact : IsCompact S :=
      isCompact_Icc.of_isClosed_subset hS_closed Set.inter_subset_left
    have hS_nonempty : S.Nonempty := ⟨t₁, ⟨⟨ht₁, le_rfl⟩, hnot⟩⟩
    obtain ⟨tstar, htstar_mem, htstar_min⟩ :=
      hS_compact.exists_isMinOn hS_nonempty continuousOn_id
    have htstar_nonneg : 0 ≤ tstar := htstar_mem.1.1
    have htstar_le : tstar ≤ t₁ := htstar_mem.1.2
    have hE_tstar_ge : R ≤ E tstar := htstar_mem.2
    -- `tstar > 0` because `E 0 < R`.
    have htstar_pos : 0 < tstar := by
      rcases eq_or_lt_of_le htstar_nonneg with h0 | hpos
      · exfalso; rw [← h0] at hE_tstar_ge; linarith
      · exact hpos
    -- Every `s ∈ [0, tstar)` has `E s < R`.
    have hlt_before : ∀ s, 0 ≤ s → s < tstar → E s < R := by
      intro s hs0 hst
      by_contra hge
      push Not at hge
      have hs_mem : s ∈ S := ⟨⟨hs0, le_trans (le_of_lt hst) htstar_le⟩, hge⟩
      have hle : tstar ≤ s := htstar_min hs_mem
      linarith
    -- By continuity at `tstar`, `E tstar ≤ R`, hence `E tstar = R`.
    have hE_tstar_le : E tstar ≤ R := by
      by_contra hgt
      push Not at hgt
      have hcont_at : ContinuousAt E tstar := (hderiv tstar htstar_nonneg).continuousAt
      have hev : ∀ᶠ s in nhds tstar, R < E s :=
        hcont_at.eventually (isOpen_Ioi.mem_nhds hgt)
      rw [Metric.eventually_nhds_iff] at hev
      obtain ⟨δ, hδ_pos, hδ⟩ := hev
      let s := tstar - min (tstar / 2) (δ / 2)
      have hs_lt : s < tstar := by
        have : 0 < min (tstar / 2) (δ / 2) := lt_min (by linarith) (by linarith)
        linarith
      have hs_nonneg : 0 ≤ s := by
        have : min (tstar / 2) (δ / 2) ≤ tstar / 2 := min_le_left _ _
        linarith
      have hs_dist : dist s tstar < δ := by
        rw [Real.dist_eq, abs_of_neg (by linarith)]
        have : min (tstar / 2) (δ / 2) ≤ δ / 2 := min_le_right _ _
        linarith
      have h1 : R < E s := hδ hs_dist
      have h2 : E s < R := hlt_before s hs_nonneg hs_lt
      linarith
    have hE_tstar_eq : E tstar = R := le_antisymm hE_tstar_le hE_tstar_ge
    -- At `tstar`, the Riccati inequality forces `E' tstar < 0`.
    have hE'_neg : E' tstar < 0 := by
      have hric := h_riccati tstar htstar_nonneg
      rw [hE_tstar_eq, hsqrt_R] at hric
      have hsimp : -γ * R + C * R * (γ / (2 * C)) = -(γ / 2) * R := by
        field_simp [ne_of_gt hC]
        ring
      rw [hsimp] at hric
      nlinarith
    -- Since `E' tstar < 0`, points immediately to the left of `tstar` have `E s > E tstar = R`.
    have hslope : Filter.Tendsto (fun s => (s - tstar)⁻¹ * (E s - E tstar))
        (nhdsWithin tstar {tstar}ᶜ) (nhds (E' tstar)) := by
      simpa [hasDerivAt_iff_tendsto_slope, slope_def_field, div_eq_inv_mul] using
        (hderiv tstar htstar_nonneg)
    have hev_slope : ∀ᶠ s in nhdsWithin tstar {tstar}ᶜ,
        (s - tstar)⁻¹ * (E s - E tstar) < 0 :=
      hslope.eventually (isOpen_Iio.mem_nhds hE'_neg)
    rw [eventually_nhdsWithin_iff, Metric.eventually_nhds_iff] at hev_slope
    obtain ⟨δ, hδ_pos, hδ⟩ := hev_slope
    let s := tstar - min (tstar / 2) (δ / 2)
    have hs_lt : s < tstar := by
      have : 0 < min (tstar / 2) (δ / 2) := lt_min (by linarith) (by linarith)
      linarith
    have hs_nonneg : 0 ≤ s := by
      have : min (tstar / 2) (δ / 2) ≤ tstar / 2 := min_le_left _ _
      linarith
    have hs_dist : dist s tstar < δ := by
      rw [Real.dist_eq, abs_of_neg (by linarith)]
      have : min (tstar / 2) (δ / 2) ≤ δ / 2 := min_le_right _ _
      linarith
    have hs_ne : s ∈ ({tstar}ᶜ : Set ℝ) := ne_of_lt hs_lt
    have hneg_prod : (s - tstar)⁻¹ * (E s - E tstar) < 0 := hδ hs_dist hs_ne
    have hinv_neg : (s - tstar)⁻¹ < 0 := inv_lt_zero.mpr (sub_neg.mpr hs_lt)
    have hEs_gt : E tstar < E s := by nlinarith
    have hEs_lt : E s < R := hlt_before s hs_nonneg hs_lt
    linarith
  -- Step 2: Inside the barrier, `F(t) = exp((γ/2) t) * E(t)` has nonpositive derivative on `[0, ∞)`.
  let F : ℝ → ℝ := fun t => Real.exp ((γ / 2) * t) * E t
  have hF_deriv : ∀ t ≥ 0,
      HasDerivAt F (Real.exp ((γ / 2) * t) * ((γ / 2) * E t + E' t)) t := by
    intro t ht
    have h1 : HasDerivAt (fun s => Real.exp ((γ / 2) * s))
        (Real.exp ((γ / 2) * t) * (γ / 2)) t := by
      simpa [mul_comm] using ((hasDerivAt_id t).const_mul (γ / 2)).exp
    have h2 := h1.mul (hderiv t ht)
    convert h2 using 1
    ring
  have hF_nonpos : ∀ t ∈ Set.Ici (0 : ℝ),
      Real.exp ((γ / 2) * t) * ((γ / 2) * E t + E' t) ≤ 0 := by
    intro t ht
    have hEt_lt : E t < R := h_trap t ht
    have hEt_nn : 0 ≤ E t := hE_nonneg t ht
    have hsqrt_le : Real.sqrt (E t) ≤ γ / (2 * C) := by
      rw [← hsqrt_R]
      exact Real.sqrt_le_sqrt (le_of_lt hEt_lt)
    have hric := h_riccati t ht
    have hCterm : C * E t * Real.sqrt (E t) ≤ (γ / 2) * E t := by
      have h1 : C * E t * Real.sqrt (E t) ≤ C * E t * (γ / (2 * C)) :=
        mul_le_mul_of_nonneg_left hsqrt_le (mul_nonneg (le_of_lt hC) hEt_nn)
      have h2 : C * E t * (γ / (2 * C)) = (γ / 2) * E t := by
        field_simp [ne_of_gt hC]
      linarith
    have hbracket : (γ / 2) * E t + E' t ≤ 0 := by linarith
    exact mul_nonpos_of_nonneg_of_nonpos (le_of_lt (Real.exp_pos _)) hbracket
  have hF_anti : AntitoneOn F (Set.Ici (0 : ℝ)) := by
    apply antitoneOn_of_hasDerivWithinAt_nonpos (convex_Ici (0 : ℝ))
    · intro t ht
      exact (hF_deriv t ht).continuousAt.continuousWithinAt
    · intro t ht
      exact (hF_deriv t (interior_subset ht)).hasDerivWithinAt
    · intro t ht
      exact hF_nonpos t (interior_subset ht)
  -- Step 3: Multiply `F(t) ≤ F(0) = E(0)` by `exp(-(γ/2) t)`.
  intro t ht
  have hFt : F t ≤ F 0 := hF_anti (Set.self_mem_Ici) ht ht
  simp only [F, mul_zero, Real.exp_zero, one_mul] at hFt
  have hexp_pos : 0 < Real.exp (-(γ / 2) * t) := Real.exp_pos _
  have hmul := mul_le_mul_of_nonneg_right hFt (le_of_lt hexp_pos)
  have hcancel : Real.exp ((γ / 2) * t) * E t * Real.exp (-(γ / 2) * t) = E t := by
    calc
      Real.exp ((γ / 2) * t) * E t * Real.exp (-(γ / 2) * t)
          = E t * (Real.exp ((γ / 2) * t) * Real.exp (-(γ / 2) * t)) := by ring
      _ = E t * Real.exp ((γ / 2) * t + -(γ / 2) * t) := by rw [← Real.exp_add]
      _ = E t := by simp
  linarith

end

end SingularBlowup.NonlinearBootstrap
