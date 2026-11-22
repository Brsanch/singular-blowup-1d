# Review of Theorem 3 (coupled (G, h) nonlinear stability), 2026-11-21

Reviewed: commit `ef0eaea` "Prove coupled (G,h) nonlinear stability with zero derivative loss (Theorem 3,
Lean 4, CH-1..4)" — `paper/stability_spectra_singular_houluo.tex`, `lean/SingularBlowup/CoupledHSector.lean`,
`validate/edge_coupled_h_nonlinear_closure.py` — and the three manuscript-only follow-up commits before it
(`afc5e06`, `3039582`, `2ed1e2a`). Line numbers are those of `ef0eaea`. At the time of the review `ef0eaea`
was a local commit on `main`, not yet pushed; this note sits on a branch cut from `2ed1e2a`.

**Review form:** one Claude session, independent of the author of `ef0eaea`, re-deriving the estimates by hand
and running the probe below. No separate reviewers ran. Not expert review. Nothing here certifies or refutes
Theorem 3.

Probe: `validate/theorem3_review_2026_11_21_probe.py` (`Q1`–`Q5`, about 15 s; `PASS` means the finding
reproduces). It does not import the reviewed gate, so it runs with or without `ef0eaea`.

## Verdict

1. The estimates of Theorem 3 hold when re-derived. No false step was found in the argument that the coupled
   energy closes.
2. The numerical verification claimed for it does not test it. Check `CH-4` of the new gate integrates a
   different system; with the manuscript's system the same run blows up, for a reason the manuscript itself
   documents (Remark `rem:w3_scope`).
3. Well-posedness for the coupled system is asserted, not proved, and three bounds the proof relies on are
   not stated anywhere. They are true.
4. Two displayed constants and one displayed equation are wrong as written. None changes the conclusion.
5. The three follow-up commits read as correct.

Severity scale: `[S1]` a reader could not tell it is wrong before relying on it; `[S2]` a gap that announces
itself; `[S3]` wrong but harmless; `[doc]` a fix that changes no mathematics.

## Findings

| ID | Severity | Where | Finding | Checked by |
|---|---|---|---|---|
| C1 | S1 (evidence) | gate lines 18, 123–124, 142–147 | `CH-4` is described as integrating "the exact coupled Gauge-(a) PDE". It treats the repo's `Q_op` as the subtracted operator `Q̃_op`, but `Q_op` keeps the edge value `H[g](1)`. The run therefore uses the Theta-balance linear operator and a velocity off by `−H[g](1)`. | probe `Q2`, reading the code |
| C2 | S1 (evidence) | tex 700–701, 621 | "verified by `edge_coupled_h_nonlinear_closure.py`, `CH-1..CH-4 PASS`". With the manuscript's right-hand side the same grid, data and time step blow up. No finite-difference run of this kind can verify the theorem. | probe `Q3`, `Q4` |
| C3 | S2 | tex 1349–1368, 1433–1436, 1455–1459 | Local well-posedness, uniqueness and continuity for the coupled system are asserted in two lines; `𝒲_h` is not defined; the bounds `‖j^{3/2} q_j‖_{ℓ²}, ‖j q_j‖_{ℓ¹} ≤ C√𝓔(G)` are not stated in any lemma; the Riccati inequality is stated pointwise although Theorem 2 was corrected to the integrated form in `3039582`. | hand derivation |
| C4 | S3 | tex 1387–1412 | The displayed `m_1`, `m_2` are the multipliers of `½Q_op`, not of `Q_op`; `p` is written with a factor `½` that contradicts its definition (line 1114); the stated `‖p_zz‖ ≤ ½·2^{1/4}√𝓔(G)` is false. The correct constant is `2^{1/4}` and the next bound still holds. | probe `Q1` |
| C5 | S3 | tex 1432, 1438–1439 | The coefficient equation for `d̃_m` omits the `p₀` term and uses the opposite sign convention to Lemma `lem:embeddings`; "once `C√𝓔(G) ≤ 2`" should be about 1.3. | hand derivation |
| C6 | doc | gate 12, 15, 34, 52, 174; tex 618–628, 1358 | Hard-coded absolute repo path; docstring says 16 test functions where 12 and 8 are used; `CH-4` uses weight 100 and far-field energies only; three overfull boxes of 54–64 pt. | reading, PDF log |

### C1, C2 — what `CH-4` integrates, and what happens with the right system

The manuscript's system (eq. `comoving_exact`, Theorem 3's proof) is

    G_τ = −2 G_z + ½ Q̃_op G + B q + p G_z + ½ (Q̃_op G) G,   Q̃_op G = Q_op G − H[g](1),   p = V_G(1) − Q_op G.

The gate (lines 142–147) sets `QtG = Q_op @ G`, `p = Vg1 − QtG − Hg1`, `dG = M_GG @ G + …`, on the stated
assumption (lines 123–124) that the repo's `Q_op` already vanishes at the edge. It does not:

- `Q_op @ G` at `z = −4` is −0.05482 against `H[g](1)` = −0.05464 (probe `Q2`), so `Q_op` is the unsubtracted
  operator. `M_GG = −2∂_z + ½Q_op` likewise.
- At `z = −8` the gate's velocity perturbation is +0.058; the manuscript's is +0.0035.
- A scratch check (not in the probe) also contradicts the comment that `M_GG = −2 D_up + ½ Q_op`: the two
  differ by up to 79.9 in an entry, with `max|M_GG|` = 159.8.

Rerun with both right-hand sides (probe `Q3`, the gate's grid `N = 800`, `z ∈ [−10, 10]`, its data and step):

| right-hand side | `G₀ = 0, q₀ = 0.2` | `G₀ = 0.15, q₀ = 0.2` |
|---|---|---|
| as written in the gate | `E_tot` ratio 3.85e-03 | `E_tot` ratio 1.10e-02 |
| manuscript's Gauge-(a) | `|G| > 1e6` at τ = 0.32 | `|G| > 1e6` at τ = 0.22 |

Mechanism (probe `Q4`). Let `L_a = M_GG − ½·1⊗h1` with `h1 @ G = H[g](1)`, the discretised `−2∂_z + ½Q̃_op`.

- `M_GG` has `max Re` = −0.001. `L_a` has exactly one eigenvalue with positive real part: +35.3 at
  `z_min = −10` and +101.8 at `z_min = −12`. The ratio 2.88 is close to `e¹`, i.e. the rate scales like
  `e^{|z_min|/2}` and diverges as the grid resolves the edge.
- The explicit eigenfunctions `G_λ` of OBS-09, sampled on the grid, satisfy `L_a G_λ ≈ λ G_λ` with relative
  residual 0.04, 0.03, 0.04 for `λ = ½, 1, 2`; under `M_GG` the residual is above 1.

So the grid version of the Gauge-(a) operator carries the growing non-analytic modes that Remark `rem:w3_scope`
describes. Theorems 2 and 3 exclude them by edge analyticity `𝓗¹(D(1,1))`, which a finite-difference grid
cannot impose. This does not contradict Theorem 3. It does mean that `CH-4` passes only because it integrates
another system, and that the sentence at lines 700–701 is unsupported.

What `CH-1`–`CH-3` do cover: `CH-1` the two multiplier bounds and a spot check of the smoothing inequality
(ratio 0.79 on 12 functions); `CH-2` the disk coefficient bounds on 8 power-law sequences, using only the
leading far-field boundary term; `CH-3` the scalar Young inequality. These are spot checks of ingredients.

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
- `edge_coupled_h_nonlinear_closure.py`: exit 0 (`CH-1..CH-4 PASS`), subject to C1 and C2.
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
