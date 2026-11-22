# Findings — hostile review of the Track B pure-step weak-edge theorem (at f5106e9)

**Date:** 2026-11-22
**Reviewed:** commit `f5106e9` — `docs/theorem_track_b_pure_step_weak_edge.md` (Theorem 1, Lemmas 1–4),
`validate/edge_track_b_proof_lemmas.py` (TB-1..4), the NM / PM / FS gates as repaired in `1a2f8ac`, and the
folding of the first review into the three weak-edge notes.
**Review form:** independent — analytic audit by one Claude session that did not write the theorem;
numerical attack by a separate Claude agent with fresh context (own real-space `Q`, own ETDRK4 solver, mutation
tests). The analytic findings were not cross-checked by a second reviewer. Not expert review.
**Scripts:** `validate/review_2026_11_22_track_b/` (the numerical agent's, unedited, not re-run by the
analytic reviewer; outputs went to a session scratch directory).
**Later commits:** `898e363`, `270b26e`, `bacc8b8` landed during the review. `270b26e` rewrites Lemma 3 with an
`s`-scaled Nishida norm, adds the `c₁` continuation and the `V` remainder, and moves the cone margin — the
same points as A1–A3 below. **Those later versions were not reviewed here.** Track A and Track C were not
reviewed.
**Does not claim:** that the theorem is false. No tracked file other than this note and the scripts folder
was changed.

Severity: **S1** a reader could not tell it is wrong before relying on it · **S2** wrong or unsupported but
announces itself · **S3** wrong but harmless · **doc** wording.

## Verdict

| Item | Verdict at f5106e9 |
|---|---|
| Lemma 1 (closedness; pole-free `Q` bound (2.1); bilinear bound (2.2); `a(𝒫)` identity (2.3) and zero-loss bound (2.4)) | holds — re-derived by hand; reproduced numerically on 16 families |
| Lemma 2 (first-line bound; rank-one split with `φ* = u/(u²+1)`; `t^{−1/4}` second-line bound (2.6)) | holds — re-derived by hand; no violation found over 16 `F` × 3 radii × `t ∈ [1e-4, 1]` |
| Lemma 3 (Nishida contraction) | the contraction estimate (2.8) is correct; the induction keeping the iterates bounded has a gap as written (A1) |
| Lemma 4 (`V = O(ε r^{5/4})`, transversal horizon) | holds — re-derived; reproduced with 30-digit quadrature and on computed solutions |
| Theorem 1 statement | three statements exceed the proof (A2, A3) |
| The constructed solution | survives: independent solver, real-space residual 8e-11 … 1.3e-7 against a nonlinear term 3e-5 … 1e-2; `c₁ = 3a(𝒫)` to 1e-5; physical-time equation with the true `V_G` to 3.5e-9 … 5.4e-8 |
| Gates TB-1..3, PM-2, PM-3, NM-3, FS-3 | do not support what they are cited for (C1–C4) |

## A. Analytic findings

**A1 [S2] Lemma 3: the bound on the iterates is not obtained by the written argument.** The norm
`sup d_a(ρ,s)^{1/2}‖·‖_{X_ρ}` has no weight in `s`. On the shrunken cone `d_{a_j}(ρ,s) ≥ (a_{j+1}−a_j)s`, so the
pointwise bound on `Ĝ^{(j+1)} − Ĝ^{(j)}` is `∝ s^{−1/2}`, not uniform as `s → 0`. The parenthetical "near
`s = 0` the integral itself is `O(ε s^{3/4})`" is not uniform near the cone boundary: the derivative-loss
term is `∫_0^s ds'/(ρ'−ρ)` with available gap `≥ Δ_j s + a_j(s−s')`, which gives `(1/a_j)·ln(1 + a_j/Δ_j)`,
i.e. `O(log j)` for the whole iterate and factorial growth if differences are bounded that way. Repair:
Nishida's weight `d/s` (the weakly singular `(s−s')^{−1/4}` term still closes under it). `270b26e` adopts an
`s`-scaled norm; not reviewed.

**A2 [S2] `c₁ = 3a(𝒫)` is asserted, not proved.** `c₁` is not defined for a general element of `X_ρ` (the
space controls `Re m ≤ 5/8`, i.e. `G = O(w^{5/8})`). It can be derived: `h = −(π/2)a(𝒫(G̃(s)))u^{−2} + O(u^{−9/4})`
for `s > 0` (continuity of `s' ↦ a(𝒫(G̃(s')))` from (2.4); the pole-subtracted forcing is `O(u^{−5/4})`
pointwise), so `Ĝ` is meromorphic on `Re m < 9/8` with a simple pole at `m = 1`, and `G = c₁w + O(w^{9/8})`.
None of that is in the f5106e9 text.

**A3 [S3] Statements stronger than the proof.** Item 2 states `V = (2/3)c₁r² + O(ε r³|ln r|)` unconditionally;
Lemma 4 proves `O(ε r^{5/4})` (enough for the conclusion), and the numerics show the true remainder for
analytic-edge profiles is `O(r³)`. "Unique" holds within the analytic scale, not in `C([0,2]; X_{ρ₁})` at one
radius. "Classically" is asserted without argument. The bound is claimed at `ρ₁ = ρ₀/2` and `s = 2`, which is
on the boundary of the cone `ρ < ρ₀ − (ρ₀/4)s`, not inside it.

**A4 [doc]** The Status line says "complete line-by-line analytic proof".

Checked and found correct, beyond the lemmas: the mild formulation is the right linear evolution for
`∂_s − ½Q` on the `a = 0` class (line-¼ conjugation); the non-emptiness example; the pointwise Hardy bounds
used for boundary terms; `∫_0^∞ x^{−5/2}(x − ln(1+x))dx = 2π/3`.

## B. What the numerics add about the theorem itself

- Gate amplitude in the theorem's units: `‖G₀‖_{X_ρ}/max|G₀|` = 120.85 / 427.3 / 3103.9 at ρ = 0.5 / 1 / π/2, so
  the gates' ε = 0.05 is `ε_thm` = 6.0 / 21.4 / 155. No gate evaluates a strip norm of the solution.
- For `max|G₀| ≤ 0.005` both line norms decay monotonically and
  `sup_s‖G(s)‖_{X_{ρ₀/2}}/‖G(0)‖_{X_{ρ₀}}` ≤ 0.60 / 0.28 / 0.08 at ρ₀ = 0.5 / 1 / π/2 — the `Cε` bound is
  reproduced at small amplitude (one data family, `h₀ = u³e^{−u}`).
- Apparent analyticity radius of the solution (spectral decay rate on the real line; biased by prefactors):
  ≈ 1.8 / 1.3 / 1.0 / 0.69 / 0.33 at `max|G₀|` = 0.005 / 0.02 / 0.05 / 0.1 / 0.2 (linear flow: π). With
  ρ₀ = π/2 the solution is outside the cone already at 0.02–0.05. This does not contradict the theorem (ε₀
  is unspecified); it bounds `ε₀(π/2)` well below the gate amplitude.
- Lemma 2 constants grow with ρ: `t^{1/4}‖S_g(t)F‖_{5/8}/|a(F)|` → 2.02 / 2.86 / 5.36 at ρ = 0.5 / 1 / π/2.
- Horizon: `s_∞` = 0.999994 … 0.970659 for `max|G₀|` = 1e-4 … 0.5, `dr/ds = −1.0000000` at the crossing.

## C. Gate findings

**C1 [S1] TB-1, TB-2, TB-3 do not test the statements they are cited for** ("lemma-by-lemma numerical
verification"). Mutation tests: TB-1 still passes with `Q` replaced by an unrelated operator (it is a
discrete summation by parts with the same `Q` on both sides); TB-2 passes with `φ* → 3φ*` and `φ* → 1/(u+1)`
(no residue or pole cancellation is computed; it runs at ρ = 0, while the constant at ρ > 0 is 1.3–3.5× larger);
TB-3 passes with `a ≠ 0` data, with `cot → tan`, and with the `1/(ζ+1)` smoothing removed, and its printed
`‖G‖_X = 125.376` is noise-contaminated (closed form 120.846). TB-4 is sound.

**C2 [S1] PM-3 / NM-3 still cannot detect a wrong `V_G`.** The new "mutation control" flips `V` only inside
the residual. With `V_G` replaced everywhere by `3V_G`, `−V_G` or `0.05·sin(7r)` both checks pass
(PM-3 residual 1.86e-7 / 1.74e-7 / 1.75e-7). A sensitive check exists and passes for the true formula
(`t5b_phys.py`: 3.5e-9 … 5.4e-8, against 2e-4 / 1.1e-4 / 3.7e-6 for the three wrong `V`).

**C3 [S1] PM-2 and FS-3 pass only on particular (domain, `N_z`) cells.** PM-2 at `N_z = 8192` on its own
`[−64,64]`: NaN. The unfiltered FS-3 solver is NaN on `[−64,64]` at 8192 / 16384, `[−80,80]` at 2048 / 16384,
`[−96,96]` at 2048 / 4096, `[−104,24]` at ≤ 8192 and `[−128,128]` at every `N_z`; `[−72,72]` is the only domain
in the scan finite at all four `N_z`. Mechanism (sponge test): the `1/u` tail of `𝒩⁻¹𝒫` wraps the periodic
seam into the undamped region; a tanh end-sponge on the forcing removes 10 of 15 NaN cells and makes domains
agree to 1e-12. The printed values themselves are right. "`[−64,64]` is clean" is not: 2e-7 bulk and 1.9e-3
edge-slope contamination by `s = 1.5`.

**C4 [S2] First-review items against the repaired gates.**

| First-review item | Status |
|---|---|
| `a = 0` true by construction, no mutation control | relabelled (PM, FS: quadrature now, still no control — with the nonlinearity replaced by `G²` PM-2 passes entirely); untouched (NM-3) |
| PM-2 `max|a|` not `dt`-converged | fixed (1.92e-8 at `N_s` = 240 / 480 / 960; domain-limited) |
| residuals are wraparound | partly fixed (domain `[−64,64]`); no `N_s` / `N_z` / domain convergence asserted; PM residual is pure `O(ds²)` |
| `N_z` convergence | not fixed — fails (C3) |
| physical-time checks blind to `V_G` | relabelled, not fixed (C2) |
| FS-3 H³ | fixed on the chosen domain (H² = 6.9862e-3, H³ = 2.5768e-2, converged on `[−72,72]`) |
| solvers NaN | untouched, broader than first reported (C3) |
| NM-1 / NM-2(c) tautologies | untouched; NM-2's `conj_err` 2.69e-9 is wraparound (4.2e-14 on `[−64,64]`) |
| FS-1 false at k = 0 | fixed by restriction to k ∈ {2, 3} |

**C5 [doc] Stale numbers.** The pure-step note's PM-2 / PM-3 lines and the cusp note's NM-3 row quote the
pre-repair runs (ε, iteration counts, residuals — including the first review's wraparound values). The FS
note says `|a_quad| < 1e-10` (gate: 1.81e-09), "exact infimum 0.8935" (gate prints 0.903016; 0.8935 is the
first review's number), "monotone-decreasing", "zero interior roughness growth" (not computed by the gate).
PM / NM / FS docstrings describe the pre-repair gates. PM-2 says "geometric convergence" while the loop hits
its cap at 29.

## D. Folding of the first review into the notes

Faithful: the cusp contraction and the finite-Sobolev closure are marked open, the `k ∈ {2,3}` restriction
is stated, the inflow-boundary framing is included, the `QG = a + Q̃G` justification replaces the wrong one.

## Not checked

Lemma 3's constants, `ε₀` and uniqueness numerically; `s ∈ (1.5, 2]`; data other than `h₀ = u³e^{−u}`; a
complex-line solve for the analyticity radius; Track A, Track C, Lean; anything after `f5106e9`; the
singular-product convention against the source paper (first review, B3).
