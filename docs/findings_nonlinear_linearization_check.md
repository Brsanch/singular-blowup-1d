# Finding: the stability operator is the correct linearization of the trusted
# nonlinear eq.(2.4) flow (decisive nonlinear cross-check; long-time march still
# edge-limited)

Date 2026-06-05. Code: `validate/nonlinear_linearization_check.py` (PASS), using the
CLM-exact-validated nonlinear engine `src/dynamic_rescaling.py`.

## What was tested and why

A full nonlinear "shape test" by long-time marching on the singular saddle fights two
walls at once (the saddle's unstable symmetry modes + the edge-resolution wall) and,
on the equispaced grid, only shows a regularization-artifact decay (inconclusive,
recorded in `findings_nonlinear_crosscheck` Part 2). The DECISIVE, feasible form of
the cross-check is to verify the stability operator IS the correct linearization of
the genuine nonlinear flow -- testing the flow's TENDENCY at finite amplitude rather
than its long-time evolution, which avoids the saddle entirely.

For a shape perturbation `(dOmega, dTheta)` supported on `X>1`:
```
T_nl := [ RHS_2.4(base + eps*dpert) - RHS_2.4(base) ] / eps        (trusted engine)
L     := analytic linearization of eq.(2.4) (assembled independently from the PDE)
```
If `||T_nl - L|| / ||L|| -> 0` (linearly in eps), the operator is the correct
linearization, catching any sign/factor/derivation error.

## Result

`rel.err` scales EXACTLY linearly with eps (the O(eps) forward-difference signature):

| eps | rel.err |
|---|---|
| 1e-2 | 4.59e-4 |
| 1e-3 | 4.59e-5 |
| 1e-4 | 4.59e-6 |

=> the stability operator is the correct linearization of the trusted eq.(2.4) RHS.
The spectral-gap / no-unstable-mode verdict rests on the right operator.

**The test earned its keep:** it first showed a flat ~47% mismatch (eps-independent),
which localized a genuine discretization inconsistency -- the base derivative
`Omega_bar_X` in the `-dU*Omega_bar_X` coupling must use the SAME upwind scheme as the
engine's transport (`d_upwind3`), not central (`d_central4`); at the sharp regularized
edge these differ by O(1). Fixing it gave the clean linear convergence above. (This is
a discretization-consistency lesson for any future implementation, not a flaw in the
derivation.)

## Independent nonlinear sign-check

Rayleigh quotient of the nonlinear tendency `<dpert, T_nl>/<dpert,dpert>` for the
shape perturbation:
- unweighted: `+0.25` -- positive, the expected NON-NORMAL transient (the unweighted
  norm is marginal/transport, cf. `findings_spectral_gap`);
- convective weighted (`gamma=0.5`): `-1.18` -- genuine decay tendency under the true
  nonlinear flow, consistent with the convective gap.

## Status

The nonlinear cross-check is complete in its decisive feasible form:
- edge/`xi` mode `lambda=1` confirmed by the exact flow (`findings_nonlinear_crosscheck`);
- the stability operator confirmed as the correct linearization of the trusted
  nonlinear RHS (this finding);
- shape decay tendency confirmed in the convective weighted norm.

Still open / edge-limited: a LONG-TIME nonlinear march watching shape decay to `Omega_bar`
needs an edge-resolving (mapped-grid) nonlinear solver with a stabilized normalization
(the equispaced march is regularization-limited). Lower priority now that the operator
is confirmed against the nonlinear flow. Not a proof.
