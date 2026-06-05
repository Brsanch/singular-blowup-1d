# Finding: resolved stability spectrum of the singular profile is stable;
# the free-boundary edge mode needs desingularization (the precise frontier)

Date 2026-06-05.  Code: `src/stability.py`, `validate/stability_spectrum.py`.
Target: Conjecture 2.4 of arXiv:2604.01868 (asymptotic stability of the singular
Hou-Luo self-similar profile).

## What was built and validated

The linearized operator about the certified singular base state (Thm 2.3), frozen
multipliers `mu=nu=0`, assembled as a 2N x 2N matrix on (omega, theta):
```
A_ww = -diag(abar) D_up - I - diag(Omega_bar_X) H
A_wt = D_c
A_tw =                      - diag(Theta_bar_X) H
A_tt = -diag(abar) D_up
```
`abar = U_bar + 2X` is the EXACT closed-form base speed (no regularization).
The singular base derivatives `Omega_bar_X ~ (X-1)^{-3/2}`, `Theta_bar_X ~
delta(X-1)` are regularized with offset `delta` / width `w`.

**Validation (no coupling):** dropping the `diag(.) H` terms gives a clean,
fully STABLE spectrum: eigenvalues `0, -1, -2, ...` (the marginal transport mode,
the `-omega` damping, etc.), `max Re = -2e-7`.  Transport (3rd-order upwind for
the outgoing `abar`), damping, and source assembly are correct.

## Result

**The leading "unstable" eigenvalue is a free-boundary regularization artifact,
NOT a genuine instability.**  Convergence study (`stability_spectrum.py`):

- delta-dependence (N=1100): `lambda_max(nontriv)` swings `-0.13 (delta=0.14) ->
  +0.86 -> +0.97 -> +6.8 -> +23.7+36i (delta=0.035)` — it DIVERGES as delta->0.
- the corresponding eigenmode localizes at the edge `x* ~ 0.997` with high-
  wavenumber fraction `~0.45` (grid-scale oscillation).

This is the `(X-1)^{-3/2}` edge blowing up like `delta^{-3/2}` under naive
regularization — exactly the free-boundary singularity not being handled.

**The RESOLVED / physical spectrum is stable.**  Filtering eigenmodes by
resolution (hi-k fraction `< 0.03`):

| N | delta | # resolved | max Re over resolved modes | where |
|---|---|---|---|---|
| 1400 | 0.06  | 340 | **-0.708** | box boundary x*=40 |
| 1900 | 0.045 | 526 | **-0.888** | box boundary x*=40 |

Every resolved mode has `Re lambda < 0`; the least-stable resolved modes are
far-field truncation modes at the box edge (x*=40), not modes of the profile
core.  No resolved unstable mode exists anywhere — and none near the profile.

## Interpretation (calibrated, not overclaimed)

- **Evidence consistent with Conjecture 2.4** in the resolved sector: the
  physical (grid-resolved) part of the frozen linearized spectrum is entirely
  stable, with no growing mode supported on the profile.
- **NOT a proof / not the full verdict.**  Two pieces remain open:
  (i) the **free-boundary edge mode** at `X=1` is unresolved by naive
      regularization (it is the delta-divergent grid artifact); a definitive
      statement about it requires the **desingularized free-boundary
      formulation** (below);
  (ii) the **scaling/symmetry modes** (`mu, nu != 0`, the normalization
      directions) are not yet included; they add specific neutral modes that
      must be projected out before reading off an asymptotic-stability verdict.

## The precise frontier: the desingularized free-boundary eigenvalue problem

The obstruction is local to the support edge `X=1`.  The correct formulation
introduces the **edge displacement** `xi(tau)` as an explicit scalar DOF and a
regular shape perturbation `g(Y)` in the edge-fixed frame `Y = X - xi(tau)`:
```
delta Omega(X) = g(Y) - xi * Omega_bar_Y(Y),   delta Theta = h(Y) - xi * Theta_bar_Y(Y).
```
The `-xi Theta_bar_Y = -xi (pi/2) delta(Y-1)` reproduces the delta-source
(the edge motion), removing it from the shape equation; the `(Y-1)^{-3/2}` term
combines with the edge-motion compatibility condition (`u(1) + mu = ...`) into a
finite `xi`-ODE.  The eigenvalue problem is then posed in the weighted space
where `g` is regular (perturbations scale like the profile `(X-1)^{-1/2}`, not
its derivative).  Building `eig` of THAT operator (edge DOF + regular shape +
closed-form U_bar) is the definitive computation; it is the analytic core that
makes Conjecture 2.4 genuinely open.

## Honest status

A genuine, validated negative-control result: the naive regularized spectrum is
dominated by an edge artifact (correctly caught by the convergence study, not
reported as an instability), while the resolved sector is stable. The remaining
work is the desingularized free-boundary eigenvalue problem, stated above at
maximum resolution.
