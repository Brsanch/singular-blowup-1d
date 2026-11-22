# Findings — the discrete `Q_op` omits the edge tail, and several gates depend on that

**Date:** 2026-11-22
**Scope:** `src/stability_edge_z.py` (`J_op`, `forcing`, `Q_op`, hence `M_GG`). Found during the hostile
review of the weak-edge notes; confirmed by the reviewing session.
**Review form:** self — the session that found the defect also wrote the option and ran the comparison.
No independent review of this note.
**Does not claim:** that any theorem in the manuscript is wrong; that either setting is "the" correct
discretization of the constrained operator.

## 1. The defect

`J_op = Cz @ diag(e^z) @ Hg_op` and the exponential average `EA` start at `z_min`. The part
`∫_{−∞}^{z_min} H(g) e^{z'} dz' ≈ a·e^{z_min}` (`a = H[g](1)`) is missing, so the discrete
`Q_op G` and `2·forcing·G` equal the continuum `Q_op G` **minus `a·e^{z_min−z}`**.

For `G = w/(2+w)²` (exact `a = −0.176777`), N = 3000, `z_min = −14`:

| | old (`edge_tail=False`) | with tail |
|---|---|---|
| `Q_op` at `z_min` | 0.00000 | −0.17673 |
| `Q_op` at `z_min+1` | −0.11153 | −0.17674 |
| `Q_op` at `z_min+3` | −0.16796 (= `a(1−e^{−3})`) | −0.1767 |
| `Hg_op` at `z_min` | −0.17674 | −0.17674 |

At fixed `z` the error vanishes as `z_min → −∞`; it is an O(1) layer within a few units of `z_min`,
independent of N. `Hg_op` is unaffected. Reading `Q_op` or `forcing` near `z_min` as the edge value `a`
returns ≈ 0 for any `G`.

## 2. The option

`build_blocks_z(..., edge_tail=False)`. `True` adds `e^{z_min}·Hg_op[0,:]` to every row of `J_op` and
the weight `e^{z_min−z}` on the edge node of `EA`. Default stays `False` (bit-identical matrices) for the
reason in section 3.

## 3. What changes with the tail on

Every `validate/` script importing the module was run twice (tail off / on), same code otherwise.
21 scripts; two timed out at 300 s in both settings (`transient_growth*.py`) and are not compared.

**Unchanged:** `spectral_gap`, `theta_munu_gap`, `edge_spectrum_continuation`, `edge_w3_conjugation_check`,
`coevolve_march`; and the headline spectrum with the inflow BC (max Re eig(M_GG) = −0.00173 / −0.00123 at
N = 1000 / 1400, zero eigenvalues with Re > 0, both settings).

**Gates that flip from PASS to FAIL:**

| Gate | What moved |
|---|---|
| `edge_eprime_probe` (Q3) | the free-edge mode is no longer at Re λ ≈ +1.7: max Re = −0.0016 |
| `edge_gauge_a_eigenmode_check` (EM-4) | `H[g](1)(T)/H[g](1)(0)` = −0.001, −0.000, −0.000 at T = 0.5, 1, 2 instead of 0.654, 0.394, 0.141 (gate expects ≈ e^{−T}) |
| `edge_nonlinear_bootstrap_closure` (NB-6) | `max‖½Q̃‖/‖G‖_{3/4}` = 2.94 against the bound 0.40 (was 0.386); NB-2 `emp_max` 0.488 → 0.900 (still under 1.40) |
| `edge_nonlinear_comoving_closure` (NL-5) | Gauge-(a) march ends at `B = 1.0075` instead of 1.1682, `‖G_z‖` 0.270 → 0.0470 instead of 0.0256 |
| `nonlinear_shape_march` | overflow; fitted rates −62.96 at every γ; the script then crashes |

**Other outputs that differ without flipping a gate:** `stability_edge_z_spectrum` (free-edge column:
+1.70 / +1.73 / +1.75 / +1.74 at `z_min` = −12 … −18 becomes −0.0010 … −0.0013), `cap_scaffold`,
`edge_layer_obstruction`, `energy_estimate_search`, `edge_coupled_resolvent_scan`, `edge_mellin_symbol`,
`edge_argprinciple_count`, `adversarial_review_2026_11_21_probe` (R1 peaks 1.079 / 4.42 → 1.089 / 4.76).
The size of those differences was not analysed.

## 4. Reading (hypotheses, not established)

- The truncated operator is the continuum one minus `a·e^{z_min−z}`: its forcing vanishes at the inflow
  node, which is compatible with the imposed `G(z_min) = 0`. The tail-inclusive operator has forcing `½a`
  at a node where `G` is pinned to 0. That mismatch is the likely reason the nonlinear marches change or
  blow up; whether the blow-up is a numerical boundary-layer effect or reflects the equation was not
  determined.
- The "spurious free-edge mode at Re λ ≈ +1.7" (manuscript Setup; `findings_edge_condition_derivation`)
  is absent when the quadrature includes the tail. On this evidence it is produced by the truncation,
  not by the absence of the inflow condition. Mechanism not derived.
- EM-4's "`H[g](1)` decays like `e^{−T}` on the grid" is a property of the truncated operator. With the
  tail, `H[g](1)` drops to ≈ 0 at once, which is what the constrained continuum semigroup does (LEM-10).
- NB-6 and the Gauge-(a) gates appear to use the truncated `Q_op` as a stand-in for `Q̃ = Q − a`. With an
  accurate `Q_op` they would need `Q̃` formed explicitly. Not checked line by line.

## 5. Not done

No gate, threshold or manuscript text was changed. The two `transient_growth` scripts were not compared.
The figures in `paper/scripts/` were not regenerated. Which setting each gate *should* use is undecided.

## 6. Reproducing

Comparison runner and outputs (session scratch, not in the repo):
`/private/tmp/claude-501/-Volumes-4TB-SD-ClaudeCode-noethersolve/24da0ae4-02c7-4541-b503-95490895dcb6/scratchpad/tailcmp/`
(`run_one.py <script> <0|1>` patches the default; `*.0.out` / `*.1.out`). To reproduce one row by hand:
`build_blocks_z(N=3000, zmin=-14, zmax=12, edge_bc=False, edge_tail=True)` and apply `Q_op` to `w/(2+w)²`.
