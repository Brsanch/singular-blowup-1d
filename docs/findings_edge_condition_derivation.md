# Finding: the edge-condition derivation for the desingularized free-boundary
# stability problem — worked out, cross-checked, and corrected vs the handoff shorthand

Date 2026-06-05. Code: `validate/probe_edge_condition.py` (load-bearing facts
checked). Target: Conjecture 2.4 of arXiv:2604.01868 (asymptotic stability of the
singular Hou–Luo self-similar profile). This is the analytic core that must be
right BEFORE the operator is assembled. It supersedes the one-paragraph sketch in
`findings_stability_spectrum.md` ("The precise frontier").

## Setup recap

Base state (Thm 2.3), `c_l=2, c_omega=-1`, `c_l+2 c_omega = 0`:
```
Omega_bar = 1{Y>1}/sqrt(Y-1),  Theta_bar = (pi/2) 1{Y>1},
U_bar(Y)  = 2 sqrt(1-Y) - 2  (Y<1),  -2 (Y>1),   abar := U_bar + c_l Y.
```
On the support side Y>1: `abar = 2(Y-1)`, so `abar(1)=0` (edge is a characteristic
stagnation point) and `abar'(1+) = c_l = 2` (hyperbolic, repelling to the right).
`U_bar` is finite & continuous (`U_bar(1)=-2`) but `U_bar_Y = -1/sqrt(1-Y) -> -inf`
as `Y->1^-` (vertical tangent). This vertical tangent is the whole story.

## 1. The moving (edge-fixed) frame is REQUIRED — the fixed-frame edge mode is a
##    linearization artifact (this is why the naive regularized spectrum failed)

Naive fixed-X-frame perturbation `delta Omega = -xi Omega_bar'` (rigid edge shift)
produces, through `u = H(delta Omega)`, the velocity-gradient `-xi H(Omega_bar') =
-xi U_bar_XX = +(xi/2)(1-Y)^{-3/2}` on `Y<1`. Integrated, the **velocity** itself
behaves like `xi (1-Y)^{-1/2} -> inf` at the edge: the edge shift appears to cost
*infinite* velocity. That divergence is **not physical** — it is Taylor breakdown
at the vertical-tangent edge.

Verified (`probe_edge_condition.py` (A)): the EXACT shifted velocity
`U_bar(Y - s)` is finite and bounded; only its linearization `-s U_bar_Y` diverges
like `(1-Y)^{-1/2}`. On the support side `Y>1`, `U_bar=-2` is flat, so a pure edge
shift adds **zero** velocity there — the coupling `u * Omega_bar_X` is finite.

=> Work in the moving frame `Y = X - eps xi(tau)`; keep `U_bar(Y-eps xi)` EXACT,
never form `-xi U_bar_XX`. In this frame the perturbation velocity is finite:
```
u tilde(Y) = U_bar(Y) + eps V(Y),   V_Y = H(g),   V(0) = -xi   (from U(0)=0),
=>  V(Y) = V_g(Y) - xi,   V_g(Y) := int_0^Y H(g) dY'.
```
This is the corrected reading of the recorded dead end "naive regularized
eigenvalues are delta-divergent edge artifacts": the artifact IS the fixed-frame
`-xi U_bar_XX` velocity divergence.

## 2. The edge ODE for xi (derived FOUR independent ways — all agree)

Ansatz in the moving frame: `Omega = Omega_bar(Y) + eps[g - ...]`, edge motion
carried by the frame; `g, h` regular (weighted `(Y-1)^{-1/2}` space). Allow
normalization multipliers `delta c_l = eps mu`, `delta c_omega = eps nu`.

**(i) Moving-frame regularity.** The O(eps) Omega-equation has the term
`[-xi_dot + V_g(Y) + (c_l-1) xi + mu] Omega_bar_Y`. `Omega_bar_Y ~ (Y-1)^{-3/2}`
is too singular for the weighted space; its edge coefficient must vanish:
```
xi_dot = (c_l - 1) xi + V_g(1) + mu.
```
**(ii) Theta delta-balance.** `Theta_bar' = (pi/2) delta(Y-1)`; collecting the
`delta(Y-1)` terms in the Theta-equation (using `abar delta' = -abar'(1) delta`,
`abar'(1+)=c_l`) gives the same `xi_dot = (c_l-1) xi + V_g(1) + mu`.
**(iii) Omega `(Y-1)^{-3/2}` balance.** Matching the `(Y-1)^{-3/2}` coefficients in
the Omega-equation gives the same ODE.
**(iv) Edge kinematics.** `X_e = 1 + eps xi`, `dX_e/dtau = (U + c_l X)|_{X_e} =
eps[c_l xi + delta U(1+)]`, with `delta U(1+) = V_g(1) - xi`:
`xi_dot = c_l xi + V_g(1) - xi = (c_l-1) xi + V_g(1)`. Same.

Key correction vs the naive `abar'(1)=c_l=2` characteristic estimate: the
self-velocity `delta U(1+) = V_g(1) - xi` contributes `-xi`, turning the
coefficient `c_l` into `c_l - 1 = 1`. **`xi_dot = (c_l-1) xi + V_g(1) + mu`.**

## 3. What "edge condition `H omega(1) -> 0`" actually means (handoff shorthand)

The handoff wrote the edge condition as `H omega(1) -> 0`. Precisely:
- `H omega = H(delta Omega) = u_X` is the velocity GRADIENT, not the velocity.
  The ODE coupling is to `delta U(1+)` (the velocity, `= V_g(1) - xi`), NOT to
  `H omega(1)`.
- The genuine "edge condition" is the **cancellation of the `(1-Y)^{-3/2}` part of
  `u_X = H omega` at the edge** (so the velocity stays finite) — the moving frame
  enforces this automatically. The leftover edge amplitude of the shape `g` is
  `B (Y-1)^{-1/2}` with `B = H(g)(1) / (2 lambda)`; THIS is where the finite value
  `H(g)(1) = H omega(1)|_support` enters. So `H omega(1)` controls the edge
  amplitude of `g`, it is not set to zero a priori.

Take-away: do not implement a literal `H omega(1) = 0` Dirichlet condition. Use
the moving frame + the `xi`-ODE; the correct edge BC is item 5 below.

## 4. The shape equations (support side Y>1, frozen mu=nu=0 core)

```
g_tau = -abar g_Y + c_omega g + h_Y - [ int_1^Y H(g) dY' ] Omega_bar_Y
h_tau = -abar h_Y
xi_dot = (c_l-1) xi + V_g(1)            [ + mu ]
```
The nonlocal term `-[int_1^Y H(g)] Omega_bar_Y` is regular: near the edge it is
`+ (1/2) H(g)(1) (Y-1)^{-1/2}` (profile scaling). `h` is a slaved transport mode.
`V_g(1)=int_0^1 H(g)` and the `int_1^Y H(g)` are nonlocal functionals of `g`
through the line-Hilbert.

## 5. Edge indicial structure + the symmetry mode (validation targets)

Local homogeneous edge solution `g ~ (Y-1)^r`, `r = -(lambda+1)/2`. Weighted-space
admissibility (`r >= -1/2`) <=> **Re(lambda) <= 0**. Hence:
- Any genuine unstable (`Re lambda > 0`) mode must have **no homogeneous edge
  singularity** — only the forced `(Y-1)^{-1/2}` part with coeff `H(g)(1)/(2 lambda)`.
  That is a **Robin-type edge BC** the operator must impose (NOT `H omega(1)=0`).
- **Predicted symmetry/normalization mode:** the pure edge DOF (`g=h=0, mu=0`) gives
  `xi_dot = (c_l-1) xi` => **`lambda = c_l - 1 = 1`**. This is consistent with the
  universal time-translation symmetry-mode eigenvalue (`= 1`) of self-similar
  blowup; together with the `mu, nu` (scaling/amplitude) directions it spans the
  symmetry subspace to be **projected out** before reading the Conjecture-2.4
  verdict. `lambda = 1` is a concrete, falsifiable check for the full operator.

## 6. Discrete consistency (checked, with one real caveat for the numerics)

`probe_edge_condition.py` (B): the discrete line-Hilbert commutes with the discrete
derivative `H(f')=(Hf)'` to 3e-15 on smooth content, and the OPERATOR commutator
`[H,D]` is machine-zero everywhere **except the two box ends** (rownorm 6.4 at
`Y=+-L`, vs 9.6e-3 at the physical edge `Y~1`, and 9e-16 after masking the
borders). Cause: the **circulant (periodic-wrap)** derivative fights the
non-periodic line-Hilbert at the box ends. Consequences:
- `H(Omega_bar')=U_bar_XX` (edge-velocity cancellation) is discretely consistent
  AT THE EDGE — the formulation is safe there.
- The box-end `O(1)` failure is exactly the source of the spurious **box modes
  (x*=40)** that currently dominate the "resolved" spectrum. Fix: use a NON-periodic
  (outflow) far-field derivative or a mapped infinite domain; never a circulant `D`.
- Near the edge use the **closed-form** `U_bar, U_bar_X, U_bar_XX`; do NOT discretely
  Hilbert-transform the singular `Omega_bar'`.

## 7. Net status and the precise next computation

The edge formulation is now fixed, cross-checked, and discretely consistent at the
edge. The operator to build (no longer a sketch):
```
unknowns (g, h on a MAPPED Y>1 grid resolving (Y-1)^{-1/2} and the 1/sqrt(Y) tail; xi scalar)
g_tau = -abar g_Y + c_omega g + h_Y - [int_1^Y H(g)] Omega_bar_Y   (+ nu Omega_bar)
h_tau = -abar h_Y                                                   (+ (mu+2nu) Theta_bar)
xi_dot = (c_l-1) xi + V_g(1)                                        (+ mu)
abar, U_bar, U_bar_X, U_bar_XX : CLOSED FORM;  H : non-periodic line-Hilbert on the map.
edge BC: drop the homogeneous (Y-1)^{-(lambda+1)/2} branch (Robin), keep forced (Y-1)^{-1/2}.
far field: outflow (no periodic wrap).
```
Validation gates before trusting any spectrum: (i) the `lambda = c_l-1 = 1` edge
symmetry mode appears and is pure-`xi`; (ii) eigenvalues converge under mesh AND
map-parameter refinement and are edge-resolved (low hi-k) with NO box modes; (iii)
cross-check the verdict against a nonlinear dynamic-rescaling DNS from a degenerate
IC (HARD STOP if exponents drift). Only after the symmetry subspace
(`xi` at `lambda=1`, plus `mu, nu`) is projected out is the sign of `max Re lambda`
the Conjecture-2.4 answer.

The remaining genuinely open analytic point handed off: whether the GLOBAL operator
(nonlocal `int H(g)` coupling, mapped infinite domain), restricted to the
weighted space and quotiented by the symmetry subspace, has any `Re lambda > 0`
mode. The local edge analysis forces such a mode to be non-singular at the edge
(encouraging for stability), but the nonlocal coupling can still in principle carry
an instability — that is the computation the corrected operator now makes
well-posed.
