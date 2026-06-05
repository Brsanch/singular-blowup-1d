# Finding: the self-similar profiles are saddles — forward marching fails,
# so stability must be computed spectrally (not by time-marching)

Date 2026-06-05.  Code: `src/dynamic_rescaling.py`, `validate/probe_rescaling.py`.

## What was built

A dynamic-rescaling solver for eq. (2.4) of arXiv:2604.01868 with the correct
structure:
- large box [-P,P], equispaced; line Hilbert (validated to 1e-16 at X=0 on
  `X exp(-X^2/2)`, exact `-sqrt(2/pi)`);
- **3rd-order upwind** transport for `(U + c_l X) d_X` (the speed `~ c_l X` is
  outgoing at both ends -> stable outflow; a periodic/Fourier method would be
  invalid because `c_l X` is not periodic);
- 4th-order central source terms;
- two normalization options: prescribe `c_l` + the X=0 constraint
  `c_omega = h0 + c_l - Theta_XX(0)/Omega_X(0)`, or fix `c_omega`.

## What was learned (decisive)

Forward time-marching from a smooth IC does **not** converge to a smooth
self-similar profile:
- with the **X=0 constraint**, `c_omega` starts at `~ +3.2` (positive growth) and
  feeds back to `+39` -> blowup;
- with **`c_omega = -1` fixed** (the known value), the field **decays to 0**.

Root cause: the system is **nonlinear through `U = H(Omega)`** (the term
`H(Omega) Omega_X` is quadratic), so a nontrivial steady profile exists only at a
**critical amplitude**.  Below it the `-Omega` damping wins (decay); above it the
stretching wins (blowup).  **The profiles are saddles.**  A naive forward march
slides off the saddle either way; only a normalization that actively projects out
the single unstable direction (or a root-finder) can sit on it.  This is exactly
why arXiv:2308.01528 / 2603.25104 use fixed-point / Newton constructions, not
time integration.

## Consequence for the deliverable (stability, Conjecture 2.4)

The right tool is **not** time-marching but the **linearized operator's
spectrum**:
1. the unstable eigenmodes *are* the saddle directions, so the spectrum directly
   answers stability (Conjecture 2.4);
2. the **singular** profile is already known in closed form and certified
   (`docs/findings_singular_steady_verified.md`), with closed-form velocity
   `U(X)=2 sqrt(1-X)-2` — so its linearized operator is directly constructible,
   no profile solve needed.

Next: assemble the linearized operator `L` about `(Omega_bar, Theta_bar)`,
`c_l=2, c_omega=-1`, on a mapped infinite grid that resolves the `(X-1)^{-1/2}`
edge and the `1/sqrt(X)` tail (the desingularization the smooth CAP framework
lacks), and compute `eig(L)`.  The dynamic-rescaling solver remains the tool for
the *nonlinear two-stage* DNS later (with a stabilized normalization).

## Honest status

The solver's components (Hilbert, upwind transport, sources) are individually
correct; the negative result (forward marching can't reach the saddle profile)
is a real, useful characterization, not a dead end — it redirects the stability
computation onto the correct (spectral) tool.
