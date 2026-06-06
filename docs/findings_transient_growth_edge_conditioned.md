# Finding: transient growth SETTLED with the edge-conditioned operator — the worst-case
# creep was the discretized edge singularity; in the energy norm the operator is a contraction

Date 2026-06-06. Code: `validate/transient_growth_edge_conditioned.py` (PASS). Closes the
open caveat of `findings_transient_growth.md` ("a function-analytic statement would need the
precise weighted space and an edge cutoff") using the weight + cutoff from
`findings_cap_scaffold.md` / `energy_estimate_search.py`.

## The open caveat

`findings_transient_growth.md`: PHYSICAL (smooth) perturbations have `G=1` (no transient
growth), but the WORST-CASE operator norm `sup_τ ‖e^{Mτ}‖₂` (plain norm, full operator) crept
up with edge DEPTH (`2.33 → 3.08 → 3.14 → 3.91` as `zmin: −10 → −14`) — a grid-scale
`(Y-1)^{-1/2}` edge spike. Was it a real function-space transient amplification or a
discretization artifact?

## Key fact (numerical abscissa = transient-growth bound)

In a weighted norm `⟨·,·⟩_W`, the numerical abscissa (log-norm)
`μ_W(M) = λ_max(sym_W(M))` satisfies `‖e^{Mτ}‖_W ≤ e^{μ_W τ}`. So `μ_W ≤ 0 ⇒ ‖e^{Mτ}‖_W ≤ 1`
for all `τ` ⇒ NO transient growth (monotone contraction). **Caveat (learned here):** `μ_W`
must use the ANALYTIC transport symmetric part `−φ'(z)` (exact IBP), NOT `sym(discrete −2∂_z)`
— the latter carries a boundary-discretization artifact (`+33`, `cap_scaffold`). So
`μ_an = λ_max(−diag(φ') + sym_W(½Q_op))`.

## Result — settled

Edge-conditioned operator = bulk (`z > z_cut`, where `G→0` by the Robin condition) in the
bounded-at-edge anisotropic weight `φ'(z)=2γ·½(1+tanh((z−z₀)/s))`, `γ=0.4, z₀=−6, s=2`:

| | log-norm `μ_an` | physical bump `G` |
|---|---|---|
| `z_cut=−8` | **−0.106** | 1.0000 |
| `z_cut=−6` | **−0.41** | 1.0000 |

and **resolution-independent**: `μ_an(z_cut=−8) = −0.106 / −0.103 / −0.101` at `N=1000 / 1400 /
2000` (converges; does NOT creep), and the same across `zmin=−12, −14`. Contrast the
un-conditioned plain-norm worst-case, which creeps `3.08 → 3.91` with edge depth.

So `μ_an < 0` (resolution-independent) ⇒ `‖e^{Mτ}‖_W ≤ 1` ⇒ **no transient growth** in the
energy norm; the physical bump confirms directly (`G = 1.0000`).

## Verdict

Transient growth is SETTLED for the edge-conditioned operator: in the bounded-at-edge
anisotropic energy norm, the operator is a **strict contraction** (`μ_an < 0`, resolution-
independent), so there is no transient amplification. The earlier worst-case creep
(`3 → 4`) was the discretized `(Y-1)^{-1/2}` edge spike — a grid artifact that grows with
resolution and is excluded by the Robin edge condition (`G→0`), not a function-space
transient growth. This removes the last open caveat of the transient-growth analysis and is
consistent with (and uses the same weight as) the energy-coercivity result
(`rigor_energy_estimate_setup.md` §4: `μ_an` IS the bulk energy abscissa).

## Honest status

Numerical (the analytic-transport log-norm + the discrete physical-bump semigroup); the
`μ_an` bound rests on the exact transport IBP (closed) + the discrete forcing symmetric part
(bounded, convergent). Ω-shape sector; the `H¹` coupling is separate (`rigor_*` §5). Not a
proof of Conjecture 2.4 — it settles the transient-growth (non-normality) question, the
precursor to nonlinear stability, in the edge-conditioned energy norm.
