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

## 5. What changes with the tail, and how far each reading holds

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

## 6. The remaining defect (not fixed)

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
