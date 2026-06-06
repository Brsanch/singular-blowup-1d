# Finding: transient growth SETTLED — physical perturbations have NONE (G=1); the
# non-normal worst-case is a grid-scale edge artifact (linear stability is robust)

Date 2026-06-05. Code: `validate/transient_growth.py` (PASS). The question beyond the
eigenvalues, and the precursor to nonlinear stability: the operator is strongly
NON-NORMAL (transport-dominated), so stable eigenvalues + a convective gap do not by
themselves bound the transient amplification `||e^{M tau}||` -- which is what can
trigger nonlinear instability. Resolved here.

## Conditioning had to be fixed first (a real trap)

The naive computation gave absurd, resolution-dependent amplification (millions, no
turnover) -- a discretization artifact: the forcing `(1/2)J/(Y-1)` is analytically
bounded (`-> (1/2)H(g)(1)`) but discretely `small/small` with `1/(Y-1)=e^{-z}` huge at
the edge. Recasting it as the bounded exponential-averaging convolution
`(1/2) int_{zmin}^z H(g) e^{z'-z} dz'` (weights `<= 1`) makes `||M||` resolution-
independent, so amplification is trustworthy.

## Result: separate physical from grid-scale

| zmin | N | G(bulk bump) | G(near-edge bump) | worst-case sup ||e^{M tau}|| |
|---|---|---|---|---|
| -10 | 800 | 1.000 | 1.000 | 2.33 |
| -12 | 1000 | 1.000 | 1.000 | 3.08 |
| -12 | 1400 | 1.000 | 1.000 | 3.14 |
| -14 | 1200 | 1.000 | 1.000 | 3.91 |

- **SMOOTH physical perturbations (bulk OR near-edge): G = 1.000 at EVERY resolution**
  -- NO transient growth, monotone decay. This is the physically relevant quantity.
- The **worst-case operator norm** creeps up, but tracks **edge DEPTH (zmin), not N**
  at fixed domain (3.08 -> 3.14 as N 1000 -> 1400 at zmin=-12). Its worst-case input
  is a **hi-k fraction = 1.00 spike at the deepest edge node** -- a grid-scale
  representation of the `(Y-1)^{-1/2}` singularity, NOT a physical perturbation.

## Verdict

The transient-growth / non-normality question is SETTLED for physical perturbations:
they have **no transient amplification (G=1)**, so the linear stability is ROBUST and
nonlinear stability is supported. The non-normal worst-case growth is confined to
grid-scale edge modes (a discretization feature of the singular edge that grows as
one resolves further into the singularity), not a physical instability. This closes
the earlier "partially open" caveat: the resolution-sensitivity was grid-scale
(hi-k=1.00), not physical.

## Status

Solid now: eigenvalue stability (two grids), convective spectral gap (specific-
perturbation decay), operator = correct linearization of the trusted nonlinear flow,
AND no physical transient growth (G=1). The frozen-sector linear stability picture is
robust and complete, including the non-normality precursor to nonlinear stability.

Not claimed: a rigorous bound on the grid-scale edge modes in a true function space
(they reflect the genuine operator singularity at the edge; a function-analytic
statement would need the precise weighted space and an edge cutoff). Still numerical;
the CAP/Lean closure remains.

**UPDATE 2026-06-06 — this caveat is now CLOSED** (`findings_transient_growth_edge_conditioned.md`,
`validate/transient_growth_edge_conditioned.py`): with the bounded-at-edge anisotropic weight +
edge cutoff (from `energy_estimate_search.py`), the numerical abscissa (analytic-transport
log-norm) `μ_an = λ_max(−diag(φ') + sym_W(½Q_op))` is NEGATIVE (`−0.10` at z_cut=−8, `−0.41` at
−6) and RESOLUTION-INDEPENDENT (converges, no creep) ⇒ `‖e^{Mτ}‖_W ≤ 1` ⇒ no transient growth
(contraction) in the energy norm. The worst-case creep `3→4` was confirmed to be the
discretized `(Y-1)^{-1/2}` edge spike, excluded by the Robin condition.
