# Findings — hostile review of the weak-edge work after f5106e9 (at b92611c), and its fixes

**Date:** 2026-11-22
**Reviewed:** commits `898e363`, `270b26e`, `bacc8b8`, `4a5bbf1`, `b92611c` — the Track A, B and C theorem notes,
`validate/edge_track_{a,b,c}_proof_lemmas.py`, the PM / NM / FS gates as repaired in `b92611c`,
`lean/SingularBlowup/WeakEdgeManifolds.lean`, and `rem:weak_edge_jump` in the paper. This is the part not covered
by `findings_session_2026_11_22_track_b_hostile_review.md` (which reviewed `f5106e9`).
**Review form:** independent — the analytic audit, the TB / TA / TCB gates and Lean were done by one Claude
session that did not write the theorems; the PM / NM / FS gates by a separate Claude agent with fresh context.
No second reviewer cross-checked either part. Not expert review.
**Post-fix review:** one further independent Claude agent reviewed the fix diff (no S1; its S2/S3 items on the TA-4 V tolerance, the |a + D| threshold and TB-5's wording are folded in).
**Fixes:** made on branch `fix/weak-edge-review-post-f5106e9` by the same session (the reviewer, not the
original author); see "Fixes" at the end. The review scripts are in `validate/review_2026_11_22_post_f5106e9/`.
**Does not claim:** that any theorem is false.

Severity: **S1** a reader could not tell it is wrong before relying on it · **S2** wrong or unsupported but
announces itself · **S3** wrong but harmless · **doc** wording.

## 1. Track B — did the post-f5106e9 commits close the earlier review?

| Item | Verdict at b92611c |
|---|---|
| A1 (Lemma 3 Nishida iteration) | **Closed.** Re-derived item 1 (constant 8√2 C_S C_P M / a), item 2 (8/3 C_* C_a M S^{3/4}), the base step (AM–GM) and the cancellation of the `s`-factors on the shrunken cone; the iterates are bounded uniformly including `s → 0⁺`. |
| A2 (c₁ = 3a(𝒫)) | **Closed, one S3.** The residue N(1)·(−π/2)·α = 3α checks. The note asserted that α(s) = a(𝒫(G̃(s))) is C¹ in s > 0; only Hölder of order > 1/4 is needed, and neither is written out. |
| A3 (statements beyond the proof) | **Mostly closed.** V = O(εr^{5/4}), the cone margin (d = ρ₀/6) and the uniqueness scope are fixed. **S3:** the "O(εr^{9/4}) for s > 0" remainder is not uniform as s → 0⁺ (Lemma 4 does not use it). |
| A4 (Status line) | **Not fixed (doc).** Still "Complete line-by-line analytic proof", listing TB-1..TB-4 only. |
| TB-1, TB-2, TB-3 | **Closed.** Mutants Q→u²Q, Q→1.001Q, φ*→3φ*, φ*→1/(u+1), a ≠ 0 data, cot→tan, drop 1/(ζ+1) all fail. **doc:** the IBP half of TB-1 holds for any Q (rel 1.6e-14 with a wrong Q); only the real-space Q comparison makes TB-1 sensitive to Q. |
| **TB-5** | **S1.** Its sup was attained at s = 0, so the ratio (0.5970) was a property of the initial data: identical under P → 0, 30P, G², a halved linear rate, or no decay. Its builtin `max` dropped NaN: a growth mutant whose whole trajectory is NaN passed. It stopped at s = 1.2, not 2. |
| PM-2 | **Partly closed, S3.** No NaN in 20 (domain, N_z) cells; the Picard loop does not converge at N_z = 16384 (difference 5e-6 on [−48,48]; on the gate's [−64,64] it stalls at a ~1e-9 round-off floor, not the iteration cap); cells differ by ~1% because the data are normalised by a grid maximum (with it fixed they agree to 1e-11); the residual accumulator is NaN-blind; data with a ≠ 0 pass (a is only measured at s = 0.6, 1.0, where it is automatically 0). |
| **PM-3, NM-3** | **S1.** A wrong V_G used consistently everywhere (3V, −V, 0.05 sin 7r) passes; the in-gate "control" only checks that the ODE and the residual use the same V. NM-3's accumulators start from `max(0.0, …)`: with the residual forced to NaN it passes and prints 0.00e+00. |
| FS-3 | **Partly closed, S2.** With the sponge, 4 of 24 cells are still NaN ([−96,96]@2048, [−128,128]@2048/4096/8192); the domain-check cell [−80,80] is the least converged finite domain; the drift accumulator is NaN-blind. |
| C5 (stale numbers) | Note numbers match the gates. **doc:** PM / NM / FS docstrings still stale (ε, tolerances, τ sets); NM's "V_G → −V_G mutation control (fixes C3)" comment mislabelled. |

## 2. Track A — the √v Θ-cusp theorem

| Item | Verdict |
|---|---|
| Lemma A1 (Φ_* = (1/2w)e^{−1/√w}, Mellin Γ(2m+2); a((½+D)f) = a(f); zero-loss IBP) | **Holds** (re-derived). |
| Lemma A2 ((2.6) integration by parts in s, algebraic ċ₀, L²_s bound for T) | **Holds.** **S3/doc:** S_{r→g}, 𝒩_r and ℛ_ρ were used but not defined. |
| **Lemma A3 (coupled contraction)** | **S1 — does not close as written.** The δŴ → δĜ block is the time integral ε∫₀ˢδŴ (no derivative loss), yet it is unbounded between the stated weights (d/s)^{1/2} and (d³/s²)^{1/2}: pointwise the bound is ε s/(γd), unbounded at the cone edge, and over the cone the ratio of norms grows like √K for data at frequency K (13, 27, 50, 88 at K = 10, 30, 100, 300; `a3_norms.py`). The claimed factor 16ε/(γρ₀) cannot hold. The Ŵ→Ŵ term ε(½+D)((QĜ)δŴ) fails the same way; rebalancing the weights still leaves a ln(s/d) loss. Theorem 2 is therefore not proved. |
| Theorem 2 item 2 (D = −r εa₀/(1 + ½εa₀s)) | **S3:** not derived in the note. Holds numerically: TA-4's edge value matches it to ≤ 4e-10 for s ∈ [0, 1.2]. |
| "classical on (0,2]" | **S3:** W is only in L²_s, so G̃_s ∈ L²_s. |
| TA-1 | **S3:** integration-by-parts identities, true for any q₀ (q₀ → 3q₀+1 passes). |
| TA-3 | **doc:** the "modified" side is the norm of one fixed function. |
| **TA-4** | **S1:** passes with V_G replaced everywhere by 3V_G, −V_G or 0.05 sin(7r); with the Term-A sign flipped; with c₀ doubled. c₀ = 0 was caught only by a division by zero. |

## 3. Track C — finite-Sobolev pure step

Lemmas C1–C4 hold (re-derived: (2.1) `‖b‖_{H²}`; the contour shift (2.2) and the sign of its residue term
−(π/2)a e^{z/2}; the common phase of the symmetrised kernel and the bounds (1+|η|), (1+|η|^{5/4}); λ_min of the
hierarchical forms 0.7320 / 0.5876 / 0.5415 for k = 2, 3, 4). **S3** in C5: "both Gelfand triples" must be one
triple on the intersection (or the two Lions solutions shown equal); the lower-norm energy identity behind
"unconditional uniqueness" for L^∞-in-time solutions is not written out. TCB gates are sound at symbol level
(no symmetrisation and a mismatched weight order both fail); **doc:** TCB-5 also passes for b ∉ H², because b is
normalised to unit H² norm, so it does not test "b barely in H²".

## 4. Lean and the paper remark

`lakelock lake build SingularBlowup.WeakEdgeManifolds` succeeds; zero `sorry`; all 8 theorems depend only on
`propext`, `Classical.choice`, `Quot.sound`. `hierarchical_u_dissipation_H3/H4` are the pointwise real quadratic
forms with the correct (∂ − ½)^m coefficients, i.e. the core of (2.8). **doc:** `nishida_contraction_bound` is
arithmetic only and its parametrisation differs from Lemma A3. **S3/doc** in `rem:weak_edge_jump`: "G ≡ 0 within
𝓗¹(D(1,1))" gives vanishing on the edge disk only, and the listed escape routes omitted the one Tracks B/C take
(non-edge-analytic profiles; c₁ = 3a(𝒫) ≠ 0 for s > 0 even from flat data).

## Not checked

The constants C_S, C_*, C_P and ε₀; uniqueness numerically; the full Lean graph (only the module target);
the R3 moment chain in the paper remark; PM-1, NM-1, NM-2, FS-1, FS-2; an N_s scan of NM-3; why the remaining
FS-3 cells blow up; whether the reference V_G formula is correct (the old gates could not tell).

## Fixes (branch `fix/weak-edge-review-post-f5106e9`)

- **TB-5** rewritten: solves to s = 2; every criterion fails on NaN; sup over (0, 2] (at this amplitude the sup sits at s = ds, so this criterion does not discriminate); the nonlinear defect must be
  nonzero and scale quadratically (d(2ε)/d(ε) = 3.976); residuals of the linear flow and of the defect equation at
  bulk points, coded apart from the solver (1.3e-5, 7.0e-5; linear also with Q_real: 1.5e-3; the linear residual uses the analytic G_lin, so it checks the generator and Q, not the solver). Grid-independent
  normalisation. Mutants P → 0, 3P, G², transport sign, linear rate/2, growth (NaN), Q → 1.01Q all fail.
- **TA-1:** adds the closed-form QG₀ vs real-space Q_real comparison (1.4e-3); q₀ → 3q₀+1 now fails.
- **TA-4:** the physical-time residual uses V_G computed independently from the reconstructed physical profile;
  mutation controls rerun the horizon ODE and reconstruction with 3V_G and −V_G everywhere (must be > 50× worse:
  142×); |V_ODE − V_indep| < 1e-8 (a 1.5% V error passed at the first tolerance, 1e-4); s-residual tolerance
  5e-8 (catches the Term-A sign flip, 6.2e-7); |a + D| < 5e-9 (catches c₀ × 1.05); ‖W(ds)‖/‖W(0)‖ ∈ (0.9, 1.1); NaN-safe accumulators. Mutants
  3V, −V, sin, V = 0, c₀×2, Term-A flip, 3·NL all fail.
- **TCB-5:** docstring states the scope above.
- **PM / NM / FS** (by the separate review agent): NaN-safe helpers replace every builtin `max`/`min`
  accumulator (also in PM-1, NM-1, NM-2, FS-1, FS-2 — FS-2's `min(1e9, nan)` let a NaN pass); grid-independent
  normalisation constants. **PM-2:** |a(G(0))| < 1e-9 with an a ≠ 0 data control that must fail; a at s > 0
  labelled automatic; supported N_z = {4096, 8192} (16384 does not converge — round-off floor, not the cap).
  **PM-3 / NM-3:** V_G in the residual from an independent adaptive quadrature on the reconstructed physical
  profile, |V_ODE − V_indep| < 1e-8, and controls rerunning ODE + reconstruction with 3V, −V, 0 (> 50×: 338×
  and 71×). **FS-3:** cells [−64,64]@4096/8192, [−48,48]@8192, [−96,96]@8192 ([−80,80] dropped); any
  non-finite value fails. Mutants now failing: wrong V in the ODE path or in the residual path, nonlinearity
  mutants in PM-3/NM-3, a ≠ 0 data in PM-2, H(0) mutants and NaN injections in NM-3, NaN cells in FS-3.
  **Still undetected:** a wrong V_G formula coded identically in both V implementations; FS-3 has no equation
  residual, so G² and P = 0 nonlinearity mutants pass it; NM-3's control margin is thin (71× vs 50×); NM-3 runs
  one (domain, N_z) cell.
- **Notes:** Track B status line and the two S3 caveats; Track A status ("Lemma A3 does not close; Theorem 2
  conjectural"), the gap stated in Lemma A3, S_{r→g} / 𝒩_r / ℛ_ρ defined, item 2 and "classical" qualified;
  Track C status and the two C5 caveats. Lean docstrings for `nishida_contraction_bound`. Paper remark wording.
- **Not fixed:** the Lemma A3 gap itself (needs new mathematics); the Hölder bound in Track B's c₁ continuation
  and the two Track C caveats are stated, not proved.
