# Review of Theorem 3 (coupled (G, h) nonlinear stability), 2026-11-21

Reviewed: commit `ef0eaea` "Prove coupled (G,h) nonlinear stability with zero derivative loss (Theorem 3,
Lean 4, CH-1..4)" — `paper/stability_spectra_singular_houluo.tex`, `lean/SingularBlowup/CoupledHSector.lean`,
`validate/edge_coupled_h_nonlinear_closure.py` — and the three manuscript-only follow-up commits before it
(`afc5e06`, `3039582`, `2ed1e2a`). Line numbers are those of `ef0eaea`. At the time of the review `ef0eaea`
was a local commit on `main`, not yet pushed; this note sits on a branch cut from `2ed1e2a`.

**Review form:** one Claude session, independent of the author of `ef0eaea`, re-deriving the estimates by hand
and running the probe below, followed by two independent Claude reviewers (see **Independent pass**). Not
expert review. Nothing here certifies or refutes Theorem 3. One finding of
the first version of this note was wrong and was caught by an automated review comment on the pull request;
see **Correction** below.

Probe: `validate/theorem3_review_2026_11_21_probe.py` (`Q1`–`Q5`, about 20 s; `PASS` means the finding
reproduces). It does not import the reviewed gate, so it runs with or without `ef0eaea`.

## Correction (same day)

The first version of this note (commit `7f4b881`) said that the manuscript's Gauge-(a) system "blows up" on
the gate's grid, that the discretised Gauge-(a) operator has a growing eigenvalue (+35 at `z_min = −10`, +102
at `z_min = −12`) carrying the OBS-09 modes, and that no finite-difference run could verify Theorem 3.
**That was wrong.** It was an artifact of this review's own probe.

- Mechanism. The grid integrates `δU` from `Y_min = 1 + e^{z_min}`, not from 1, so the repo's `Q_op` satisfies
  `Q_op G = H[g](1)(1 − e^{z_min − z})` near the edge (probe `Q2`: −0.03457 against −0.03456 at `z = −9`). The
  first probe subtracted the full `H[g](1)` from this truncated operator. That leaves a forcing
  `−½H[g](1)e^{z_min − z}` concentrated at `z_min`, which feeds back through the weight `e^{−z/2}` in `H[g](1)`
  with gain `~e^{|z_min|/2}`.
- Test. With the tail restored, `Q_exact = Q_op + (e^{z_min}/(Y−1))⊗h1` — the discretisation that
  `edge_nonlinear_comoving_closure.py` and check `EM-4` of the eigenmode probe already use — the largest real
  part is −0.0011 at `z_min = −10` and −0.0010 at `z_min = −12`. Without the tail it is +35.3, +101.8, +282.4
  at `z_min = −10, −12, −14`, the predicted `e^{|z_min|/2}` scaling (probe `Q4`).
- Consequence. With the tail restored the manuscript's system decays on the gate's grid (probe `Q3`). The grid
  does not carry the OBS-09 modes: under the tail-corrected flow `H[g](1)` of `G_1` decays like `e^{−T}`
  (0.654, 0.393, 0.141 at `T = ½, 1, 2`), where the continuum mode grows like `e^{+T}`. This is what `EM-4`
  had already recorded; the first version of this note contradicted it.
- Credit. The error was pointed out by the Cursor Bugbot review comment on pull request 4, which named the
  missing `Q_exact` tail and `EM-4`. It was confirmed here by the experiment above before anything was changed.

Findings C1 and C2 below are restated accordingly. C3–C6 are unaffected.

## Independent pass (same day, after the correction)

Two independent Claude reviewers with fresh context reviewed `ef0eaea` without seeing this note: one checked
the statement and proof of Theorem 3 line by line; one checked the well-posedness claims and the gate and Lean
evidence. The coordinating session read their reports and, for the two claims new to it, their outputs; it did
not re-run their scripts. Their scripts are not in the repo.

Agreement. Both found no error that changes the Riccati inequality, the decay rates or the trapping of `B`,
and both independently reported C1, C3, C4 and C5. The second reviewer also met the blow-up of a naive
correction of `CH-4`, diagnosed it as a grid artifact, and obtained the same corrected values as probe `Q3`
(ratios 4.6e-3 and 1.5e-2; `B(3)` = 1.026 and 1.082).

What they add:

- **`CH-4` is uninformative even when corrected.** With every nonlocal and nonlinear term deleted (pure
  transport) the check still passes (ratios 0.084 and 0.107 against the threshold 0.829); so does the system
  with the sign of `½Q̃G` flipped (0.0046, 0.060). The initial data are Gaussians in `z`, not analytic at the
  edge, so they lie outside `𝒳 × 𝒳_h`; `c_l` and `c_ω` are printed, not asserted.
- **Why the gate's system damps `H[g](1)` so fast.** As coded, `H[g](1)` decays at rate
  `1 + (2/3π)e^{−z_min/2}` (predicted −12.59, −32.49, −86.61; measured −12.58, −32.48, −86.58 at
  `z_min = −8, −10, −12`), where the manuscript's identity gives −1.
- **`M_GG` uses central differences**, not the upwind form the gate's comment states (Frobenius distance 0 to
  the central form, 2765 to the upwind form).
- **Step 1:** the displayed `m_1` also has the wrong sign (the true value at 0 is `+1/π`).
- **Step 3:** the sharp joint inflow bound is `1.669·Z_{h,disk}` and the dissipation left is `1.50·D_{h,disk}`
  (this note's cruder count gives 1.3; the manuscript implies 2). `‖j²q_j‖_{ℓ²} ≤ √Z_disk + √(Y_far/π)`.
- **Statement:** uniqueness is claimed with no class; `C_tot` does not depend on `B₀`; lines 633–637 say
  Theorems 2 and 3 are "proved analytically in full" while the coupled well-posedness is one sentence.
- **Lean header** cites a "Remark B.1 (`rem:h_sector`)" that no longer exists.

Severity. Both reviewers rate the evidence sentence at lines 700–701 `[S1]`, and the first rates the false
`p_zz` bound `[S1]` as a displayed bound with no effect on the conclusion. C1 and C2 below carry the
reviewers' rating; C4 is listed as S3 with that dissent recorded here.

## Verdict

1. The estimates of Theorem 3 hold when re-derived. No false step was found in the argument that the coupled
   energy closes.
2. Check `CH-4` of the new gate does not integrate the system it names. A run of the manuscript's system on
   the same grid also decays, so the qualitative outcome stands; the gate's printed values belong to another
   system.
3. Well-posedness for the coupled system is asserted, not proved, and three bounds the proof relies on are
   not stated anywhere. They are true.
4. Two displayed constants and one displayed equation are wrong as written. None changes the conclusion.
5. The three follow-up commits read as correct.

Severity scale: `[S1]` a reader could not tell it is wrong before relying on it; `[S2]` a gap that announces
itself; `[S3]` wrong but harmless; `[doc]` a fix that changes no mathematics.

## Findings

| ID | Severity | Where | Finding | Checked by |
|---|---|---|---|---|
| C1 | S1 (evidence) | gate lines 18, 123–124, 142–147 | `CH-4` is described as integrating "the exact coupled Gauge-(a) PDE". It treats the repo's `Q_op` as the subtracted operator `Q̃_op`. `Q_op` is the unsubtracted operator cut at the grid edge, so the run uses the Theta-balance linear operator and a velocity off by about `−H[g](1)`. The manuscript's system, discretised with the tail, also decays, with different numbers. | probe `Q2`, `Q3`; reading the code |
| C2 | S1 (evidence) | tex 700–701, 621 | "verified by `edge_coupled_h_nonlinear_closure.py`, `CH-1..CH-4 PASS`" is not supported. `CH-4` passes for pure transport and for a sign-flipped operator, uses data outside the theorem's space, one grid and `T = 3`; and the grid does not carry the non-analytic Gauge-(a) modes, so no grid run tests the edge-analytic hypothesis. | probe `Q4`; independent reviewer's runs |
| C3 | S2 | tex 1349–1368, 1433–1436, 1455–1459 | Local well-posedness, uniqueness and continuity for the coupled system are asserted in two lines; `𝒲_h` is not defined; the bounds `‖j^{3/2} q_j‖_{ℓ²}, ‖j q_j‖_{ℓ¹} ≤ C√𝓔(G)` are not stated in any lemma; the Riccati inequality is stated pointwise although Theorem 2 was corrected to the integrated form in `3039582`. | hand derivation |
| C4 | S3 | tex 1387–1412 | The displayed `m_1`, `m_2` are the multipliers of `½Q_op`, not of `Q_op`; `p` is written with a factor `½` that contradicts its definition (line 1114); the stated `‖p_zz‖ ≤ ½·2^{1/4}√𝓔(G)` is false. The correct constant is `2^{1/4}` and the next bound still holds. | probe `Q1` |
| C5 | S3 | tex 1432, 1438–1439 | The coefficient equation for `d̃_m` omits the `p₀` term and uses the opposite sign convention to Lemma `lem:embeddings`; "once `C√𝓔(G) ≤ 2`" should be about 1.3. | hand derivation |
| C6 | doc | gate 12, 15, 34, 52, 174; tex 618–628, 1358 | Hard-coded absolute repo path; docstring says 16 test functions where 12 and 8 are used; `CH-4` uses weight 100 and far-field energies only; three overfull boxes of 54–64 pt. | reading, PDF log |

### C1 — what `CH-4` integrates

The manuscript's system (eq. `comoving_exact`, Theorem 3's proof) is

    G_τ = −2 G_z + ½ Q̃_op G + B q + p G_z + ½ (Q̃_op G) G,   Q̃_op G = Q_op G − H[g](1),   p = V_G(1) − Q_op G.

The gate (lines 142–147) sets `QtG = Q_op @ G`, `p = Vg1 − QtG − Hg1`, `dG = M_GG @ G + …`, on the stated
assumption (lines 123–124) that the repo's `Q_op` already vanishes at the edge like `O(e^z)`. It does not:
`Q_op G = H[g](1)(1 − e^{z_min − z})` there (probe `Q2`), and `M_GG = −2∂_z + ½Q_op` likewise. So the gate's
linear and quadratic terms use `Q_op` where the system has `Q̃_op`, and at `z = −8` its velocity perturbation
is +0.058 where the manuscript's is +0.011. A scratch check (not in the probe) also contradicts the comment
that `M_GG = −2 D_up + ½ Q_op`: the two differ by up to 79.9 in an entry, with `max|M_GG|` = 159.8.

Rerun on the gate's grid (`N = 800`, `z ∈ [−10, 10]`), data and time step (probe `Q3`):

| right-hand side | data | `E_tot` ratio at `T = 3` | `B` | `c_l` |
|---|---|---|---|---|
| as written in the gate | `G₀ = 0, q₀ = 0.2` | 3.85e-03 | 1.0006 | 2.0003 |
| as written in the gate | `G₀ = 0.15, q₀ = 0.2` | 1.10e-02 | 1.0017 | 2.0003 |
| manuscript's Gauge-(a), tail restored | `G₀ = 0, q₀ = 0.2` | 4.57e-03 | 1.0260 | 2.0027 |
| manuscript's Gauge-(a), tail restored | `G₀ = 0.15, q₀ = 0.2` | 1.52e-02 | 1.0818 | 2.0054 |

All four are far below the gate's threshold `e^{−T/16}` = 0.829. The decay that `CH-4` reports therefore also
holds for the system it claims to integrate; the values it prints (`B`, `c_l`, `c_ω`, the ratios) are those of
a different system. The repair is to use `Q_exact` and `Q̃ = Q_exact − H[g](1)` as
`edge_nonlinear_comoving_closure.py` does.

### C2 — what the gate can and cannot show

`CH-1` checks the two multiplier bounds and spot-checks the smoothing inequality (ratio 0.79 on 12 functions).
`CH-2` checks the disk coefficient bounds on 8 power-law sequences, using only the leading far-field boundary
term. `CH-3` checks the scalar Young inequality. `CH-4` is one grid, two initial data, one time horizon.
These are spot checks of ingredients and one illustration, not a verification of a nonlinear stability
theorem.

A grid run also cannot test the hypothesis that makes the theorem true. Remark `rem:w3_scope` shows the
Gauge-(a) linearisation has growing solutions in every finite-order edge class; the theorem excludes them by
edge analyticity. The truncated grid does not carry those modes at all (probe `Q4`, and `EM-4`): the
tail-corrected discrete operator is stable, and `H[g](1)` of the explicit mode `G_1` decays under its flow. So
decay on the grid says nothing about whether analyticity is needed, in either direction.

### C3 — what the proof needs and does not state

The proof of Lemma `lem:embeddings`(3) already contains the decomposition

    q_j = −(1/(π j(j+1))) Σ_k (k−½) c_k/(k−j−½) + G(0)/(2π j(j+1)(j+½)) − (1/(π(j+1)(j+½))) ∫₀^∞ e^{−(j+½)z} G_z dz.

The first term is controlled by the Hilbert matrix (norm π); the last by the Carleman–Hilbert inequality
`Σ_j |∫₀^∞ e^{−(j+¼)z} φ dz|² ≤ π‖φ‖²` with `φ = e^{−z/4}G_z`. Hence `‖j² q_j‖_{ℓ²} ≤ C√𝓔(G)`, which gives
both bounds used in Step 3 and the boundedness of the zeroth-order coefficient `1 − ½p − p_z` on `𝒳_h` that
any well-posedness argument needs. The disk and far-field parts are each `±G(0)/(π j²)` to leading order, so
neither `j^{3/2} q_j^{disk}` nor `j^{3/2} q_j^{far}` is in `ℓ²` unless `G(0) = 0`: the cancellation is
necessary and should be stated. The lemma as stated gives only `Σ√j |q_j|`.

Also missing: a definition of `𝒲_h`; the truncation/mollification argument for the `q` energy (the `q_z`
equation has a coefficient `p_zz ∈ L²` only); the integrated form of the Riccati inequality in the theorem
statement.

### C4 — Step 1 constants

By the manuscript's own formula `𝓜[Q_op G](s) = 2c(s+½)Ĝ(s+½)` (line 1249), the multiplier of
`∂_z²(Q_op G)` acting on `(G_z)^(s+½)` is `2s²c(s+½)/(s+½)`, twice the displayed `m_2`. Since
`sup|2m_2| = 1`, the displayed inequality `‖∂_z²(Q_op G)‖_{L²(dz)} ≤ ‖G_z‖_{L²_{1/2}}` is correct and sharp
(probe `Q1`: ratios 0.9989, 0.9995, 0.9773). With `p = V_G(1) − Q_op G` this gives
`‖p_zz‖ ≤ ‖G_z‖_{L²_{1/2}} ≤ 2^{1/4}√𝓔`. The stated `½·2^{1/4}√𝓔` = 0.595√𝓔 fails: `Q1` measures 0.953√𝓔
and 0.978√𝓔 for narrow bumps supported just right of `z = 0`. The commutator bound that follows needs
`2√(3/2) ≈ 2.45 ≤ 3`, so `3·2^{1/4}√𝓔 (X_h + Y_h)` stands.

### C5 — Step 3 details

With `p(1+v) = p₀ − Σ_{j≥1} q_j v^j` (Lemma `lem:embeddings`), the `q` equation gives

    ∂_τ d̃_m + (2m+1) d̃_m = (m+½) p₀ d̃_m − (m+½) Σ_{j+k=m, j≥1} q_j d̃_k.

The extra term is bounded by `3|p₀| D_{h,disk}`, the same form as the others. The inflow terms cost
`2√2·Li₂(¼) + 2√2/3 ≈ 0.757 + 0.943 = 1.70` times `Z_{h,disk}`, so after reserving rate 3 only
`1.3·D_{h,disk}` of the dissipation `4D_{h,disk} + 2Z_{h,disk}` is left for the nonlinear terms.

## Checked and found correct in Theorem 3

- The `q` equation `∂_τq + (2−p)q_z + (1 − ½p − p_z)q = 0` and the forcing `+B q` in the `G` equation.
- Linear damping `−3(X_h + Y_h)` on `[z₀, ∞)` and `−4D_{h,disk} − 2Z_{h,disk}` on the disk; rate 3 for `𝓔_h`
  after the inflow.
- The Mellin identity on the pole line `Re s = 0`: `∂_z` removes the constant `H[g](1)` that distinguishes
  `Q_op` from `Q̃_op`, and `s·c(s+½)` is regular at `s = 0`.
- `‖G_z‖²_{L²_{1/2}(ℝ)} ≤ (2/3)Z_disk + √2·Y_far ≤ √2·𝓔`.
- The Sobolev bound `‖e^{−z/4}q‖²_∞ ≤ (3/2)(X_h + Y_h)` and the commutator estimate.
- `‖q‖²_𝒳 ≤ 20𝓔_h`; the forcing term `2B⟨G, q⟩_𝒳`; the Young inequality with `μ_h = 1000B₀²` (also the Lean
  lemma); the rates 1/16, 1/32, 1/64; the trapping of `B(τ)` in `[½B₀, 2B₀]`.
- The class `h(1+v) = Σ d_k v^{k+½}` is preserved by the flow.

## The follow-up commits `afc5e06`, `3039582`, `2ed1e2a`

Read in full; no error found.

- W2: the separate half-line memberships (`L²_σ ⊂ L²_{σ−½}` on `(−∞, 0]`, the reverse on `[0, ∞)`) are right
  and repair the earlier claim that `e^{−σz}G ∈ H¹(ℝ)`.
- Theorem 2 and Proposition C.1 now state the Riccati inequality in integrated form for the solution.
- The weak-edge remark now limits "no nontrivial solution" to the linearised flow.
- The new sentence attributing the weighted transient of `M_GG` to the edge value `H[g](1)` is confirmed
  (probe `Q5`): peaks 1.167 at τ = 0.0125 (γ = 0.5) and 12.727 at τ = 0.0075 (γ = 0.75) with an independent
  propagator, and `G(z = −10, τ = 0.0025)` = −4.20e-04 against `½H[g](1)τ` = −4.55e-04. The earlier wording
  "O(e^z) edge tail" was wrong.
- The numbers added for `W2-C` (Parseval error 1.42e-16, rate 18.8498) match the gate's output.
- One mismatch: the text gives `τ ∈ [0, 6]` for `transient_growth.py` (lines 47, 250). That script samples
  `[0, 5]` in steps of 0.2; `transient_growth_edge_conditioned.py` is the one that uses `[0, 6]`.

## Builds and gates

- `lakelock lake build` on `ef0eaea`: success (8257 jobs); 0 `sorry`, 0 `axiom`; vacuity lint 0 findings. The
  five lemmas in `CoupledHSector.lean` are real-number arithmetic, as the manuscript says.
- `edge_coupled_h_nonlinear_closure.py`: exit 0 (`CH-1..CH-4 PASS`), subject to C1.
- PDF from `ef0eaea`: 0 errors, 0 undefined references, three overfull boxes (60 pt and 64 pt in the gate-list
  paragraph at lines 618–628, 54 pt in the display at line 1358).
- The reproducibility gate list and both earlier probes, run on `917a496`: 14 of 14 exit 0. No gate file
  changed between `917a496` and `2ed1e2a`.

## Not checked

- A written-out well-posedness proof for the coupled system.
- The manuscript against arXiv:2604.01868.
- Any change after `ef0eaea`.
- An independent delta audit of `bc9620b` (the fold of the first adversarial review) was started with three
  reviewers on 2026-11-21 and stopped before any of them reported.

## Reproduce

    reslock run --cpu 2 --purpose "theorem 3 review probe" -- python3 validate/theorem3_review_2026_11_21_probe.py
