# Findings — the discrete edge average omitted its tail; fix, and what two reviews established

**Date:** 2026-11-22
**Scope:** `src/stability_edge_z.py` (`J_op`, `forcing`, `Q_op`, hence `M_GG`), the gates that use them, three
manuscript passages.
**Review form:** independent — two Claude reviewers with fresh context, one per version of the fix. The first
broke version 1 (section 2). The second verified the version-2 kernel and gates, and found that the
conclusions first drawn from version 2 overreached (sections 5–6); this note states only what survived.
Reviewer scripts: `validate/review_2026_11_22_qop_tail/{v1,v2}/` (not re-run by the implementer).
Not expert review.
**Does not claim:** anything about the Gauge-(a) operator, the edge problem across `X = 1`, Conjecture 2.4,
or a converged value of any edge-layer eigen/sup quantity.

## 1. The defect

`J(Y) = ∫_1^Y H(g) dY'` and the bounded form `Q_op G = J/(Y−1)` were accumulated from `z_min`. The piece
`∫_0^{v_min} H(g) dy` (`v_min = e^{z_min}`) was missing, so the discrete `Q_op G` and `2·forcing·G` equalled
the continuum value minus `(tail)·e^{−z}`: the first row of `Q_op` is identically 0 and the error is O(1)
within a few units of `z_min`, independent of N. For `G = w/(2+w)²` (exact edge value `a = −0.176777`,
N = 3000, `z_min = −14`): `Q_op` = 0.00000 / −0.11153 / −0.16796 at `z_min` / `z_min+1` / `z_min+3`.
Four gates already compensated by hand (`Q_exact = Q_op + (e^{z_min}/(Y−1)) ⊗ h1`), so the truncation was
known to them but undocumented in `src/`.

## 2. Version 1 was wrong (record)

It approximated the tail by `H(g)(z_min)·e^{z_min}` (the Hilbert row at the first node). Reviewer 1: that row
is harmonic near the pinned node, O(1) wrong on edge-scale structure and divergent like `ln(1/dz)`; the
"stiff eigenvalue" it produced was N-dependent; a fourth hand-tail gate (`edge_coupled_h_nonlinear_closure`,
indirect importer) was double-counting while printing PASS. Superseded.

## 3. Version 2 (default since 2026-11-22; Bryan: "Make default")

For `g = v^{−1/2}G` supported on `v ≥ v_min`,

```
∫_0^{v_min} H[g](1+y) dy = −(1/π) ∫_{z_min}^{z_max} e^{z'/2} · (−ln(1 − e^{z_min−z'})) · G(z') dz'      (exact)
```

`edge_tail_row(z)` integrates this kernel against the grid's hat functions (Gauss–Legendre per cell, adaptive
quadrature on the first, log-singular, cell). `build_blocks_z(..., edge_tail=True)` adds it to `J_op` (hence
`forcing`, `M_GG`) and `Q_op`.

Verified by reviewer 2 (own derivation and scripts):

| Item | Result |
|---|---|
| kernel formula, sign, Hilbert convention | correct; independent double quadrature agrees to 1e-11 – 1e-13 on four test functions |
| `edge_tail_row` weights | O(dz²); exact on constants and linears to 5e-9; first/last weights match mpmath to 1e-12 |
| `edge_tail=False` | bit-identical to the pre-change code on every returned array, 18 grid/BC combinations, four β |
| `forcing` vs `½Q_op`; β ≠ 0 conjugation | consistent (relative ≤ 6e-9; 9e-12) |
| `Q_op G` accuracy, bulk data (`w/(2+w)²`, Gaussian at z = 1) | relative 0.9e-3 – 7.4e-3, decreasing with N |
| headline spectrum with the inflow BC | unchanged |
| stiff edge eigenvalue of `M_GG` | a complex pair, −18.2 ± 23.0i / −141.3 ± 153i at `z_min` = −10 / −14, N-independent; it is the truncated operator's existing pair (−12.7 ± 18i / −95 ± 111i) scaled by ≈ 1.48, not a new eigenvalue |
| explicit time steps | `dt·ρ ≤ 0.97` at `dt = 0.25·dz` on every marching grid except (1800, −16, 14) of `nonlinear_shape_march`, where it is 2.36 and that script's `2/ρ` cap binds; no instability found |

## 4. Gate repairs (all verified byte-identical to the pre-change output by reviewer 2, except as noted)

| Gate | Cause | Repair |
|---|---|---|
| `edge_nonlinear_comoving_closure`, `edge_nonlinear_bootstrap_closure`, `edge_gauge_a_eigenmode_check` (EM-4), `edge_coupled_h_nonlinear_closure` (CH-4) | they add a tail by hand; the default double-counted it | `edge_tail=False` (`coevolve_moving_frame.build` gained an `edge_tail` passthrough) |
| `edge_eprime_probe` Q3 | asserts the free-edge `+1.7` mode of the truncated average | Q3 pins `edge_tail=False`; Q3b records that with the tail this matrix has no unstable free-edge mode at its grid |
| `nonlinear_shape_march` | version-1 stiffness | `dt ≤ 2/ρ(M_GG)`; printed output identical with the tail on and off |
| `transient_growth_edge_conditioned` | > 25 min, not re-baselined | pinned to `edge_tail=False`; identity of its output not verified by a run |

The hand tail in those gates uses the edge functional `a` instead of the exact kernel; the two differ on
structure within O(1) of `z_min` (hand/exact = 0.77 on a bump one unit above `z_min`, 0.54 on a jump).
Six scripts, including `edge_coupled_h_nonlinear_closure.py`, hard-code the main checkout as `ROOT`; a copy in
another worktree therefore imports this checkout's `src/` (reviewer 2 hit this in `.claude/worktrees/`).

## 5. What changes with the tail, and how far each reading holds  *(state before the Hilbert-diagonal repair; see section 11 for the current default)*

Sweep: every `validate/` script importing the module directly or through `coevolve_moving_frame`, plus
`review_2026_11_21_w3/disc.py`, tail off and on. All exit 0 in both settings.

**Unchanged output:** `coevolve_march`, `coevolve_paper_match`, `coevolve_stress`, `edge_spectrum_continuation`,
`edge_w3_conjugation_check`, `nonlinear_shape_march`, `stability_edge_z_spectrum` (BC, whole-spectrum maximum),
`theta_munu_gap`, the pinned gates.

### 5a. Free-edge spectrum — discretization-defined; no conclusion

This matrix, rightmost eigenvalue over `|Im| < 20`, N = 1600:

| `z_min` | free edge, truncated | free edge, with tail | BC, truncated | BC, with tail |
|---|---|---|---|---|
| −6 | +1.375 | −0.066 | −0.809 | −0.814 |
| −10 | +1.638 | −0.406 | −0.958 | −0.977 |
| −14 | +1.733 | −0.722 | −0.952 | −1.106 |
| −18 | +1.745 | −0.865 | −0.837 | −1.022 |

Reviewer 2: with the tail there is no `Re > 0` eigenvalue for `z_min ≤ −4` on any grid tried (N = 30–2400,
`z_max` = 3–30), but there is for `z_min ≥ −3` (+0.224 at `z_min = 0`); the `|Im| < 20` value is a window cut
on a continuous arc and moves with `z_max` and N; and an independent, convergent collocation of `Q_op`
(`v2/exactmat.py`) has a real unstable eigenvalue at the edge node with a free edge (+16.4 / +21.7 / +27.6
at N = 400 / 800 / 1600). A free inflow boundary has no continuum limit, so the free-edge spectrum is +1.7,
nothing, or +log N depending on the discretization. What stands: the truncated average is where the
long-quoted `+1.7` comes from; with the inflow BC both versions have `max Re < 0`.

### 5b. Field of values of the nonlocal forcing — the sign for γ ≤ ¼ survives; the magnitudes do not

`cap_scaffold`, weighted bound `C_force(γ)` of the forcing in `L²(e^{−2γz})` (its grid: N = 1200, `z_min = −12`):

| γ | truncated | with tail | independent collocation, same grid (reviewer 2) |
|---|---|---|---|
| 0.10 | +35.15 | +1.2e-4 | −5.0e-6 |
| 0.20 | +31.94 | +2.7e-4 | −5.0e-6 |
| 0.25 | — | +3.7e-4 | −5.0e-6 |
| 0.30 | +29.75 | +0.789 | +0.164 |
| 0.40 | +29.64 | +15.2 | +9.81 |
| 0.50 | +56.32 | +94.3 | +83.9 |

Established (reviewer 2, independent matrix, N = 400–3200, `z_min` = −6 … −16, `z_max` = 7 and 12, BC on and
off): the largest eigenvalue of the symmetric part of the weighted forcing is ≤ 0 through γ = 0.26 and > 0
from γ = 0.27, and for γ > ¼ it diverges with depth like `e^{−z_min/2}`. That is the sign pattern of manuscript
Proposition `prop:energy_form` (threshold ¼; the grids resolve it only to "between 0.25 and 0.27").

Not established, and not to be quoted:
- the repo-matrix values at γ ≥ 0.3 (4.8× / 1.55× / 1.12× the independent values, and N-dependent);
- "≈ 0" at γ ≤ 0.25 for the repo matrix: its small positive values shrink with N at fixed depth
  (γ = 0.2, `z_min = −10`: 2.2e-4 / 1.5e-4 / 5.8e-5 / 1.9e-5 at N = 400 … 3200) but grow with depth at fixed N
  (γ = 0.25, N = 1400: 2.9e-4 / 2.7e-3 / 2.4e-2 / 0.21 at `z_min` = −12 / −16 / −20 / −24);
- `energy_estimate_search`'s `A = −0.500` at γ = 0.25 means only `C_force(0.25) < 5e-4` on that one grid;
- `edge_layer_obstruction`'s `A_full ≈ 0` is not a γ test: it uses the ramp weight
  `φ′ = 2γ·½(1+tanh((z+6)/2))`, γ = 0.4 (exponent 0 at the edge), and its sign flips with N
  (−0.025 … +0.074 for N = 700 … 4000 at `z_min = −10`).

What can be said about the old record: the `+30 … +56` field-of-values bounds and the `+40`-type edge
abscissas (OBS-01, OBS-03, `findings_cap_scaffold.md`, `edge_layer_estimate.md`) were dominated by the
truncated edge average at every γ, including the ramp weights. They are not evidence about the operator.

## 6. The remaining defect (not fixed)  *(REPAIRED later the same day — section 11)*

`hilbert_support_z` sums over `[Y_0, Y_max]` but puts the principal value of a constant over `[1, Y_max]` on
the diagonal, so row i assumes `g ≡ g_i` on `[1, Y_0]`, while the tail assumes `g = 0` there. On structure
within O(1) of `z_min`, `Hg_op` — hence `Q_op`, with or without the tail — is 14–37% wrong, independent of N
(bump at `z_min+1`: 0.156; jump at `z_min`: 0.35–0.37; `v2/t2.py`, `t2b.py`). It predates this change.
Reviewer 2's proposed repair: remove `g_i·ln(v_i/(v_i−v_0))/π` from the diagonal for i ≥ 1; in its test this
brings the bump error to ≈ 1% (converging) and reproduces the independent collocation's `C_force` to under
1%. Not applied here: it changes every downstream number again and needs its own review.

## 7. Other changed outputs (tail off → on; not analysed)

- `edge_eprime_probe`: Q2 residuals `r` = 0.586 … 32.6 → 5.04 … 276.7; Q4 `alpha_edge_bc` 3.0285 → 3.1502
  (`findings_session_2026_07_22_eprime_residual_floor.md` quotes the old KEY line).
- `adversarial_review_2026_11_21_probe` R1 peaks at γ = 0.75: 4.42 → 6.01 and 12.73 → 17.07; the γ = 0.5 peak
  1.167 → 1.163.
- `review_2026_11_21_w3/disc.py`: every line changes; with the tail the value equals `e^{4γτ}` at τ = 1 for
  every γ and grid (hypothesis: reversed group velocity of central-difference sawtooth modes; not inspected).
  The "Discrete operator" paragraph of `findings_session_2026_11_21_w3_conjugation.md` §4 describes the
  truncated operator.
- `coevolve_projected` ("nu only" final norm), `edge_coupled_resolvent_scan` (deep Θ→G gain 0.9892 → 0.9994),
  `edge_mellin_symbol`, `edge_argprinciple_count`, `spectral_gap` (4th digit).

## 8. Manuscript

Three passages changed: Setup and Sensitivity item 1 now say the `+1.4 … +1.7` mode comes from the truncated
edge average and that the free-edge spectrum is discretization-defined; the Limitations energy item says the
earlier large positive bounds came from the truncation and reports only the sign pattern of 5b. Not changed:
the "2.33–3.91" worst-case sentence (from `transient_growth.py`, which builds its own truncated operator and
still prints those numbers). Figures not regenerated.

## 9. Not done

The Hilbert-diagonal repair (section 6). `transient_growth*.py` on the exact tail. Desk documents quoting
truncated-operator numbers (CORR notes and lowered status only on OBS-01, OBS-03). Figures. The narrative
text printed by `cap_scaffold`, `edge_layer_obstruction`, `energy_estimate_search` (each now prints a NOTE
that it describes the truncated operator). The hard-coded `ROOT` in six scripts.

## 10. Reproducing

`python3 validate/stability_edge_z_spectrum.py` (5a), `python3 validate/cap_scaffold.py` (5b, repo column),
`validate/review_2026_11_22_qop_tail/v2/exactmat.py` + `t5*.py` (independent collocation; the scripts
reference the reviewer's scratch layout and may need path edits), `python3 validate/edge_eprime_probe.py`.

## 11. Hilbert-diagonal repair (later 2026-11-22; Bryan: "Yes, repair it")

**Review form:** independent — a third Claude reviewer with fresh context (scripts in
`validate/review_2026_11_22_qop_tail/v3/`, not re-run by the implementer). No further review round was run
on the wording changes that fold its findings.

**Change.** `hilbert_support_z(..., grid_support=True)` takes the principal value of the constant over
`[Y_0, Y_max]` for rows ≥ 1, so the Hilbert quadrature and the tail both model `g = 0` below the grid.
`build_blocks_z(hilbert_grid_support=None)` follows `edge_tail`, so every gate pinned to `edge_tail=False`
is unaffected (reviewer: 1056 array comparisons bit-identical to the pre-change code).

**Verified by the reviewer**

| Item | Result |
|---|---|
| derivation of the diagonal; half-weight first node | correct; `Hz @ 1` on rows ≥ 1 equals the exact transform of the indicator of `[Y_0, Y_max]` to 1.8e-15 |
| accuracy on grid-supported data (bumps at `z_min+1`, `z_min+½`; jump at `z_min`; zero-edge-functional pair), N = 400 → 3200 | `Q_op`: 1–4% → 0.1–0.5% (jump: 6–10% → 1–2%), first order; historical diagonal 15–43%, flat in N |
| headline spectrum with the inflow BC | `max Re` and `#Re>0` unchanged; rightmost with `|Im| < 20`: −0.929 / −1.083 / −1.167 (truncated / tail only / current), independent collocation −1.159 |
| stiff edge pair | −182.9 ± 80.6i (N = 1400), −183.4 ± 81.2i (N = 2800) at `z_min = −14`; collocation −183.4 ± 81.6i |
| explicit steps | `max|R(dt·λ)| ≤ 0.99995` on all 11 marching grids; `dt·ρ` = 0.50–0.93 (1.43 on one N = 300 grid; 2.26 on (1800, −16, 14), where the shape-march cap binds) |
| scripts | 23/23 exit 0; ten byte-identical to the previous default, incl. the pinned gates, `nonlinear_shape_march`, `spectral_gap` |

**Where the new default is WORSE (S1, not fixed — a modelling trade-off).** For data with `g(Y_0) ≠ 0`,
i.e. `G ~ v^{1/2}` at the edge (e.g. `G = v^{1/2}/(1+v)^{3/2}`), `Q_op` is 2–12% wrong and `Hg_op` 20–100%
wrong, independent of N and decreasing only like `1/|z_min|`; the historical diagonal with the tail is
≤ 0.3% there. Mechanism: the removed diagonal term is exactly the contribution of `g ≡ g(Y_0)` on `[1, Y_0]`.
For such data use `edge_tail=True, hilbert_grid_support=False`. The marches (Gaussian at z = 1) are
unaffected, and no current gate shows a consequence.

**Field of values of the weighted forcing, current default** (largest eigenvalue of the symmetric part,
`L²(e^{−2γz})`, inflow BC, `z_max = 12`):

| γ | `z_min = −10`, N = 800 / 1600 / 3200 | `z_min = −12` | `z_min = −14` |
|---|---|---|---|
| ≤ 0.25 | 7.0e-7 / 4.8e-7 / 3.1e-7 | 5.9e-7 | 5.3e-7 |
| 0.27 | 7.3e-5 / 8.4e-5 / 9.0e-5 | 1.2e-3 | 4.5e-3 |
| 0.30 | 0.0593 / 0.0595 / 0.0596 | 0.165 | 0.452 |
| 0.40 | 3.56 / 3.56 / 3.56 | 9.85 | 27.1 |

Reviewer's qualifications, all folded into the manuscript sentence:
- the γ ≤ ¼ residual is carried by the last nodes at `z_max` (scales like `e^{−z_max/2}`; 7.5e-6 at
  `z_max = 7`); with the top three nodes removed the value is negative and tends to zero
  (−4.9e-8 … −1.2e-11 for N = 400 … 6400). It no longer grows with depth (5e-7–7e-7 from `z_min` = −10 to −28);
- the values for γ > ¼ depend on `z_max` (γ = 0.3, `z_min = −10`: 0.046 / 0.0595 / 0.0614 at `z_max` = 7 / 12 / 20);
- the first positive γ moves toward ¼ with depth (0.27 at `z_min = −10`, 0.264–0.266 at −14, 0.26 at −22);
- an independent second-order collocation agrees to ≤ 0.5% for γ ≥ 0.3, to 1–17% at γ = 0.27, and has the
  opposite (negative) sign for γ ≤ 0.26; it uses the same "g = 0 below the grid" model, so it checks the
  quadrature, not that modelling choice.

**Free edge.** With the current default the free-edge matrix has an edge-node mode whose growth rate increases
with N (4.4 ± 22i, 10.8 ± 30i, 18.2 ± 40i, 26.8 ± 53i for N = 400 … 3200 at `z_min = −12`); consistent row-0
choices give a larger real mode. Together with section 5a this confirms the free-edge spectrum is
discretization-defined. `edge_eprime_probe` Q3b is now an informational print of three discretizations
(+1.7066, −0.0016 + 105.2i, +13.09 + 33.1i) with no assertion; a pinned regression value was not added.

**Printed numbers that moved again** (previous default → current): `cap_scaffold` `C_force` 1.2e-4 / 2.7e-4 /
0.789 / 15.2 / 94.3 → 0.0000 / 0.0000 / 0.165 / 9.85 / 84.0 at γ = 0.1 … 0.5 — and it now flags γ = 0.3 as
"DISSIPATIVE; best", which is a one-grid statement (`C_force(0.3)` = 3.36 at `z_min = −18`);
`edge_eprime_probe` `r_min` 5.04 → 6.38, `r(−14)` 276.7 → 346.7, `alpha_edge_bc` 3.1502 → 3.4473 (beyond that
gate's own 2% dual-run tolerance for its KEY line); `edge_mellin_symbol` part (2) coefficients by under 1%;
`edge_coupled_resolvent_scan`, `coevolve_*`, `disc.py` in the 3rd–6th digit.

**Pre-existing, noted by the reviewer, not addressed:** the top end of the Hilbert quadrature is inconsistent by
half a cell (O(dz)); `Y = 1 + e^z` loses the node spacing to round-off for `z_min ≲ −28`.

**Superseded statements in this note:** section 5a's "with the tail there is no Re > 0 eigenvalue for
`z_min ≤ −4`" and the "with tail" columns of 5a/5b describe the previous default (tail with the historical
diagonal, still available as `hilbert_grid_support=False`); section 6 is repaired; section 9's first item is done.
