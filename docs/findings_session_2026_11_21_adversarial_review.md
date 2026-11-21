# Findings — adversarial review of the manuscript (2026-11-21)

**Date:** 2026-11-21
**Reviewed:** `paper/stability_spectra_singular_houluo.tex` and its gates at commit `0dbf217`. Line numbers below refer to that commit.
**Probe:** `validate/adversarial_review_2026_11_21_probe.py` (R1–R3, exit 0, about 15 s) reproduces the three numerically checkable findings.
**Review form:** four independent Claude reviewers with fresh context and no access to earlier review notes, one per lens (linear theory; nonlinear estimates; well-posedness; evidence versus claims), each on a frozen snapshot of `0dbf217`. The coordinating Claude session (which had itself edited parts of the manuscript earlier that day) re-ran E1, E2 and N1's base identity and re-read the code or text for E3, E5, E8, part of E4, N4 and the P/L items. Everything else rests on the reviewers' own checks. This is not expert review.
**Does not claim:** that either theorem is false · that either theorem is proved · anything about the PDE beyond the specific statements below.

**Status against `main`.** One commit landed during the review (`2508e65`). Where a finding was visibly addressed there, the last column says so; items not marked were not re-checked against `2508e65`.

---

## 0. Verdict in five lines

1. Inside its stated scope, no reviewer found a false statement in Theorem 1 or Theorem 2.
2. Two numerical-evidence claims are wrong as written (E1, E2).
3. Theorem 2 is a theorem about the support-interior reduction. No nontrivial solution of it is shown to satisfy the original system across the edge `X = 1`, and the "codimension-1 manifold" on which it would is not invariant (N1).
4. Five proof steps are unsupported as written (P1, P2, L1, L2, N2). The reviewers believe all five are repairable.
5. Several gates test less than the text says (E3–E7).

## 1. Findings

Severity: **S1** a reader could not tell it is wrong before relying on it · **S2** wrong or unsupported but announces itself · **S3** wrong but harmless · **doc** wording.

### S1

| ID | Where | Claim as written | What the evidence shows | At `2508e65` |
|---|---|---|---|---|
| E1 | l.116–117, 273, 286–287, 723; `nonlinear_shape_march.py:105,127,150–151` | "peak ‖·‖_w/‖·‖_0 = 1.000" for the finite-amplitude march | The gate computes the peak in the **unweighted** norm (γ = 0), at τ = 0.5, 1, …, 3 only, floored at 1. On a fine time grid the linear operator `M_GG` (inflow BC on, the gate's bump) gives, in `L²(e^{−2γz}dz)`: γ = 0.5 → 1.08, 1.17, 1.18 and γ = 0.75 → 4.42, 12.73, 20.90 for `z_min = −10, −14, −16`. The peak is N-independent (12.73 at N = 1400, 12.72 at N = 2000) and sits at τ ≈ 0.003–0.05. With nonlinearity on, the reviewer measured 1.166 / 1.167 / 1.173 at γ = 0.5 for A = 1e-3 / 0.5 / 4. `spectral_gap.py` fits rates on τ ∈ [0.5, 2.5], which skips the transient. Probe R1. | still present (l.117, 286–287, 704, 728) |
| E2 | abstract l.45–46; l.243–248; l.666–668; `transient_growth_edge_conditioned.py` | "zero physical transient growth in the edge-conditioned energy norm"; only a "≈3.3% short-time amplification" is quoted; "numerical abscissa strictly negative and resolution-independent" | The gate's own column `G_mat = sup_τ‖exp(M_b τ)‖_W` is 111.94 at N = 1000, `z_min = −12`, `z_cut = −8`, attained at the last sampled time τ = 6 and still rising (234 at τ = 8). The discrete log-norm is +33.4 / +47.0 / +67.3 for N = 1000 / 1400 / 2000. The reviewer located the top eigenvector at the outflow end `z_max` (one-sided stencil), not at `z_cut` as l.243–244 says; that attribution was not re-checked. Smooth bumps do give 1.0000, and `μ_an` reproduces. Probe R2. | partly addressed: the Sensitivity section and the gate now say "sampled matrix norm ~112 … not a discrete matrix contraction"; the abstract sentence is unchanged, and "~112" is a sample at τ = 6, not a supremum |

### S2

| ID | Where | Claim as written | What the evidence shows | At `2508e65` |
|---|---|---|---|---|
| N1 | l.549–551, 641–643 (and l.143) | The weak edge relation `B = A²` holds on "the codimension-1 manifold `H[g](1) = 0`"; only the case outside it is open | `{H[g](1) = 0}` is not invariant. For the linear Gauge-(a) flow `d/dτ H[g](1) = −H[g](1) − c₁/3` (`c₁ = G_v(0)`), because `∫ e^{−z/2}·½Q̃_op G dz = (π/3)c₁`. More generally, with `F_m = f.p.∫ w^{−m−3/2} G dw`, `Ḟ_m = −(2m+1)F_m + (π/(2m+3))c_{m+1}` (reviewer: verified numerically for m = 0, 1, 2). So `H[g](1,τ) ≡ 0` forces every `c_m ≡ 0` and every far moment to vanish, i.e. `G ≡ 0` at linear order. Consequence: no nontrivial solution of Theorem 2 is shown to solve (eq:dyn) across `X = 1`, and the limit `(Ω̄, B_∞Θ̄)` with `B_∞ ≠ 1` is not a weak steady state by the Setup's own criterion (l.169–173). The remark's edge balances and the defect `(π/2)(B−1)δ` re-derive correctly. Probe R3 (the m = 0 identity). | still present (l.143, 550, 643) |
| P1 | l.509–513, 1251–1263, 1283–1285, 1370–1376, 1418 | "applying the identical three-block estimates … gives `d/dτ 𝓔 ≤ …`" for the iterates and for the solution | Raised independently by two reviewers. Step 1 gives a mild solution in `C([0,T]; 𝒳)`. The estimates differentiate `𝓔` and use `D_disk = Σ m³|c_m|²` and `G_zz`, neither finite for a general element of `𝒳` (`𝒳` allows `D_disk(0) = ∞` and jumps of `G_z` at `v > 1`). No regularisation or integrated form is given. Believed repairable: truncated sums `Z_N` close on the disk because both convolutions are lower-triangular; Friedrichs commutator on `[z₀,∞)` with `∂_z p ∈ L^∞`; pass to the limit by weak lower semicontinuity. Not written. | — |
| P2 | l.1414–1415 | `G ∈ C([0,T₀]; 𝒳)` "follows from weak-* continuity together with continuity of the energy via the differential identity of Lemma B.4" | Lemma B.4 contains an inequality, not an identity, and it is not established for the limit. An upper bound on `d𝓔/dτ` plus weak continuity gives right-continuity only. Time reversal is unavailable (the disk flow is dissipative). Needs an integrated energy identity (which needs `∫D_disk < ∞` from P1). | — |
| L1 | l.757–762 | "By Lemma B.1(1), `(λ−M_σ)G = F` on `Re s = σ ∈ (0,½)` transforms to the ladder equation" | The citation does not apply. (a) Lemma B.1(1) is about `Q̃_op`; `M_σ` uses `Q_op` (real-space quadrature: `Q_op ↔ c(w)` on `Re w ∈ (−½,½)`, `Q̃_op ↔ c(w)` on `Re w ∈ (½,3/2)`). (b) B.1(1) covers `Re s ∈ (0,¼]` and `G ∈ L²_{3/4,1/4}`. (c) `Ĝ(s+½)` is not a Mellin transform in general: for wave packets `F̂ = e^{−(τ−τ₀)²/2}e^{iξτ}` at σ = 0.3, λ = 0.4 the reviewer found `|Ĝ(σ+½+iτ₀)|/|F̂(τ₀)| = 20.1, 40.0, 80.0` at τ₀ = 10, 20, 40, which indicates `R_slit: L²_σ → L²_{σ+½}` is unbounded (medium confidence). Not load-bearing: the right-inverse paragraph and Lemma A.4 each derive what they need. | — |
| L2 | l.858–862 | "showing `G ∈ D(M_σ)`" by shifting the contour of `c(w)Ĝ(w)` from `Re w = σ+½` to `σ` | No bound on `Ĝ` in the open strip or on `Re s = σ+½` is given; Lemmas A.1–A.3 bound only the `Re s = σ` trace. Needed: `c·Ĝ ∈ H²(strip)`, i.e. `(λ+2s)R_slit` bounded on `L²`. Two repairs: rerun A.1–A.3 with the extra factor `(λ+2s)/(λ+2s′)` (log-derivative ≤ 1/δ), or prove it for `F̂ ∈ C_c^∞` and extend by closedness. The conclusion is believed true. | — |
| N2 | l.1291–1309 (Remark h-sector) | The augmented energy `𝓔 + μ_h𝓔_h` satisfies the Riccati inequality | The linear part checks out (factor 20; rate −3 including inflow terms 0.347 and 0.472 absorbed by the disk rate −6). The nonlinear part needs bounds Lemma B.3 does not supply: `j^{3/2}`-weighted control of `q_j` for `(∂_z p)F`, and `∂_z²(Q G) ∈ L²_{−1/4}` in the far field. `B₀` is undefined. Asserted, not proved. | — |
| E3 | l.594–596; `edge_nonlinear_comoving_closure.py:283–304` | Gate check NL-6 "weak-edge jump balance `B = A²`" | The check integrates a fixed Lorentzian against one test bump and confirms the result is π/2. Its second assertion is algebraically the first. No solution, `A(τ)`, `B(τ)`, `H[g](1)` or PDE enters. More broadly, "gates verify all nonlinear identities and constants" overstates: `C_nl` is never computed, no gate evaluates `𝓔` along a trajectory, and NB-5 integrates a scalar ODE with an assumed constant. | — |
| E4 | l.117–119, 258–262, 289–298; `coevolve_projected.py`, `coevolve_moving_frame.py` | Co-evolving march: "ξ ≡ 0 preserved"; relaxation to (2, −1) | (a) `ξ` never enters `dG` or `dh`, and with `xi_mu = −Vg1` the code is exactly `ξ′ = (c_l−1)ξ`; `ξ(0) = 0` stays 0 by construction, and `ξ(0) = 1e-8` grows to 5.6e-7 by τ = 4. (b) The gates print "h-force only: BLOWS UP", "both (+feed c_l): BLOWS UP", "forcing (μ,ν): BLOWS UP (saddle runaway)"; `coevolve_paper_match.py` says Theta-balance was "chosen for numerical stability". None of this is in the paper. (c) The march uses `c_ω = −c_l/2`, which is neither Gauge (a) of Theorem 2 nor eq. (3.2) of the source; by l.539 that normalisation gives `Ȧ/A = ½A·H[g](1) ≠ 0`, incompatible with the imposed `G(z_min) = 0`. (d) "verified against the raw PDE to 2.8e-15" is one `(G,h)` pair with both sides built from the same discrete operators. | — |
| E5 | abstract l.61–70; l.132–135, 633–636 | Abstract states the nonlinear theorem without the `h ≡ 0` restriction | Theorem 2 is stated "in the invariant Θ-sector `h ≡ 0`" (l.506); `h ≢ 0` exists only as Remark h-sector (N2). | still present |
| E6 | l.447–452; `edge_w1_route_f_closure.py:112–113, 233` | "F1–F5 PASS at σ = 0.4, δ = 1.1" as verification of W1 | That anchor is `Re λ = +0.30`; F5 runs at σ = ¼, λ = 0.3. No W1 or W2 gate samples a point with `−2σ < Re λ < 0`. W2-A/B test one function at one `s`. The quoted "max 9.3587" is a 28-sample maximum (the supremum tends to 3π ≈ 9.4248, still below 15.41). | — |
| E7 | l.44–45, 104, 209–210, 662; `src/stability_edge.py`, `verify_singular_steady.py` | "two independent discretizations", both "with the inflow condition" | `stability_edge.py` contains no inflow-condition code path. Its rightmost eigenvalue is −0.94 (box modes), not ≈ −1e-3. Its discrete `H(Ω̄)` is a sawtooth of +11.6 / +4.0 / +2.5 at near-edge nodes where the exact value is 0 (reviewer's run). `verify_singular_steady.py` prints numeric `U(1) = −1.686` and a Hilbert error of 0.149. | — |

### S3 and doc

| ID | Where | Note |
|---|---|---|
| N4 | l.1186–1187, 1156 | `(4/3)·∫_{−∞}^{−ln 2} e^{3z/2}dz = 2√2/9`, not `4√2/9`. The stated bound is still valid (twice the true value). |
| L3 | l.980–982 | "`e^{−(σ−½)z}G₊ ∈ H¹([0,∞))` so `lim e^{−sz}G = 0` for all `Re s ≥ σ−½`" is a false inference when `−2σ < Re λ ≤ 1−2σ` (counterexample `G = χ(z)e^{−λz/2}`). Harmless. `2508e65` rewrote this passage; not re-reviewed. |
| L4 | l.998–1000 | "`Ĝ` has `L²_loc` boundary values on `Re s = σ−½`" fails when `Re λ = 1−2σ` (the pole is on the line); `P = Ĝ/X` is unaffected. |
| L5 | l.401–405, 436 | The Remark says the constraint `Ĝ(½) = 0` is essential in Lemma A.4. For the theorem's operator (`T_c` on `Re s = σ`) Lemma A.4 never uses it: `P(0) = 0` comes from the `Γ(s)` pole alone. The sentence is true only of the Gauge-(a) operator. l.436 writes `½Q̃_op G = e^{−z/2}T_c G`, contradicting l.358. Theorem and Remark are consistent: `(λ−M_σ)G_λ = Γ(1+λ)·1 ∉ L²_σ`. |
| L6 | l.450 | "`‖X‖₂‖K_slit‖₂`: 0.3441 → 0.0405"; accurate quadrature gives 0.33886 and 0.04107. |
| L7 | l.176–182, 357 | `M_GG = −2∂_z + ½Q_op` is not the linearisation at fixed `(c_l, c_ω)`; that one carries an extra `½e^{−z}V_G(1)`. It is the linearisation in the gauge `c_l = 2 − V_G(1)`, `c_l + 2c_ω = 0` (`B` frozen, `A` free). "Frozen-exponent" is a misnomer. Noted by two reviewers and by the coordinator. |
| P3 | l.1338–1340, 1365–1367 | `a^{(n)} ∈ [3/2, 5/2]` on `z > 0` needs `0.32√R₀ ≤ ½`, not stated. Pullback by `Φ_ℝ` is bounded `𝒳 → 𝒳`, not on the far component alone. The contraction factor is about 0.62, not ½. Conclusions unaffected. |
| P4 | l.509, 1323, 1413 | "Solution in `C([0,T]; 𝒳)`" is never defined. Uniqueness is proved among solutions with `𝓔 ≤ R₀`; the theorem states it without that restriction. |
| E8 | l.665, 724, 704, 206 | Stale `(2.000006, −1.000003)` and "NL-1..5" (both fixed in `2508e65`); "{1200..2200}×{−14,−16}×{7,10}" is five runs, not a product. |
| E9 | `NonlinearBootstrap.lean:350–357` | `riccati_bootstrap_trapping_decay` assumes the inequality for all `t ≥ 0` and a globally differentiable `E`; the manuscript proves it only on `{𝓔 < R}`, so global existence is an input to the Lean theorem, not an output. |
| misc | l.374, 534, 544, 585, 850, 878, 936 | `D` is evaluated on the slit (should be `D₊`); "Rankine–Hugoniot" is a misnomer; "classically on `X > 1`" holds only in the strong `L²` sense for `H¹` data; "tr = 2.14799" is a truncation of 2.14803; Stirling is uniform only on compact sets; the "smooth bounded kernel" has a jump at τ′ = 0; `c₀` names two different constants; `T_c` is never defined. |

## 2. The structural point behind N1, E1 and E4(c)

Three reviewers arrived at pieces of one picture, which matches `OBS-09`:

- The weak form of the Ω-equation across `X = 1` requires `B = A²` (edge amplitude of Θ equals the square of the edge amplitude of Ω).
- **Gauge (a)** pins `A ≡ 1` and lets `B` drift: `Ḃ/B = −H[g](1)`. Theorem 2 lives here. Its energy estimate is sound, but `B ≠ 1` leaves a point-mass defect at the edge.
- **The Theta-balance gauge** (`c_l + 2c_ω = 0`, used by the numerics and equal to the operator of Theorem 1) pins `B` and lets `A` drift: `2Ȧ/A = A·H[g](1)`. Then `G` does not stay zero at the edge, i.e. the solution leaves `L²_σ` for σ > 0. On a truncated grid this shows up as the weighted-norm transient of E1, which grows with edge depth.
- In either gauge, staying inside the ansatz `Ω = (X−1)^{−1/2}(1+G)·1_{X>1}`, `Θ = B(π/2 + h)·1_{X>1}` with `B = A²` needs `H[g](1) ≡ 0`, and by N1 that is not preserved by the flow.

So neither theorem, as it stands, is a stability statement for solutions of the original system. What would be needed is what the manuscript's own Remark names: an edge boundary layer in Θ or a dynamic edge exponent. This paragraph is the reviewers' and coordinator's reading, not a proved statement.

## 3. Checked and found correct

- **Derivation of (eq:comoving_exact)** from (eq:dyn) under Gauge (a): all formulas and signs (two reviewers, by hand).
- **Lemma B.1:** multiplier, range `Re s ∈ (0,¼]`, constants 2/5 and 1/2.
- **Lemma B.2:** disk Hilbert expansion; coefficient formulas against an independent contour computation (1e-8); `‖𝓗_{1/2}‖ = π`; the `b_m` bound and series.
- **The three linear energy inequalities,** tested directly from the Hilbert transform on a 16-function span including kinked and growing far fields: no violation; worst `𝓔̇/𝓔 = −0.755` against the claimed −1/8.
- **Lemma B.3:** every constant, with `1/√π` sharp for the sum norm; pointwise bounds; the `ℓ¹_{1/2}` bound.
- **Lemma B.4:** every boundary term at `z₀`; all four terms of `∂_z𝒩`; `|𝒩_disk| ≤ 3(|p₀| + 2Σ√j|q_j|)D_disk` for both nonlinear terms; the Lyapunov rates; the Riccati step.
- **Slit-plane resolvent:** orientation, signs, `κ_split`, `κ₁`, the `D(1)` subtraction; `R_slit` reproduces two exact solutions to 2e-23 at four parameter points. Closed form found: `R_slit(λ)X(·;μ) = [X(·;μ) − Γ(1+μ)/Γ(1+λ)·X(·;λ)]/(λ−μ)`.
- **Lemmas A.1–A.3:** margins, exponents, `L(σ,δ)`, both `L²` bounds; Schur sums ≤ 0.45 for `|ω| ≤ 300`.
- **Lemma A.4:** continuation, pole orders including λ = 0, Morera periodicity, the Parseval factor ½, the exponent `p`, the exclusion `8π|k| < 6π`.
- **Remark on eigenmodes:** contour range, both far-field asymptotics, `H[g_λ](1) = −2Γ(1+λ)`, `r₁, r₂, r₃`.
- **Appendix C:** `𝒳` is a Hilbert space with no derivative loss; `S` bounded `𝒳 → 𝒳`; the disk flow is contractive on `H¹(D)`; every term of the weak-norm contraction, with no derivative loss; `T₀` depends only on `R₀`.
- **Gates:** all quoted numbers other than those in E1, E2, L6, E8 reproduce.

## 4. Not checked by anyone

- Whether (eq:dyn), Theorem 2.3 and the singular-product convention match arXiv:2604.01868 line by line.
- A written repair of P1, P2, L2 or N2.
- The full `transient_growth.py` run and `G_mat` on the other three grids.
- `lake build` and `#print axioms` (not run in this review).
- Any time-dependent simulation of the coupled disk and far-field system in `𝒳`.
- The changes in `2508e65` (Lemma A.4 for `0 < δ ≤ 1`, gate conclusion, Section 9–10 summaries).

## 5. Reproducing

- R1 (E1), R2 (E2), R3 (N1, m = 0): `python3 validate/adversarial_review_2026_11_21_probe.py`.
- E3: read `validate/edge_nonlinear_comoving_closure.py` lines 283–304.
- E4(a): `validate/coevolve_moving_frame.py` line 97, `dxi = (c_l − 1)·xi + Vg1 + xi_mu` with `xi_mu = −Vg1`.
- E6: `−2σ + δ = −0.8 + 1.1 = +0.3`.
- The reviewers' own scripts were written to a session scratch folder and are not kept in the repo. The parameters above are enough to regenerate each number.
