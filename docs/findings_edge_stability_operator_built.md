# Finding: desingularized free-boundary stability operator BUILT — frozen-sector
# spectrum has no genuine unstable mode (numerical evidence for Conjecture 2.4)

Date 2026-06-05. Code: `src/stability_edge.py`, `validate/stability_edge_spectrum.py`
(PASS), building on the verified edge-condition derivation
(`docs/findings_edge_condition_derivation.md`, `validate/probe_edge_condition.py`).
Target: Conjecture 2.4 of arXiv:2604.01868 (asymptotic stability of the singular
Hou–Luo self-similar profile).

## What was built

The desingularized linear-stability operator about the certified singular base
state (Thm 2.3), `c_l=2, c_omega=-1`, frozen multipliers `mu=nu=0`. Three structural
moves make it well-posed where the naive operator produced a delta-divergent edge
artifact (recorded dead end):

1. **Moving (edge-fixed) frame** — the fixed-X-frame edge shift `-xi Omega_bar'`
   has a SPURIOUS infinite edge velocity (Taylor breakdown at the vertical-tangent
   edge); that was the artifact. The moving frame keeps the perturbation velocity
   finite.
2. **`xi` drops out of the shape equation** — substituting the edge ODE
   `xi_dot=(c_l-1)xi+V_g(1)+mu` turns the singular edge term into
   `-[int_1^Y H(g)] Omega_bar_Y`, a pure functional of `g`. The operator is
   **block-triangular**; the full spectrum is `eig(M) U {c_l-1}`. The pure-edge
   (`xi`) symmetry mode at `lambda = c_l-1 = 1` is appended ANALYTICALLY, so the
   edge singularity is never put on the grid.
3. **Weighted variable** `G=(Y-1)^{1/2} g` — because `c_omega=-1` exactly cancels
   the weight's contribution, the `G`-equation has no zeroth-order term and a
   BOUNDED forcing:
   `G_tau = -abar G_Y + (Y-1)^{1/2} h_Y + (1/2) J(Y)/(Y-1)`, `J=int_1^Y H(g)`,
   with `(1/2)J/(Y-1) -> (1/2) H(g)(1)` finite at the edge.

Closed-form `U_bar, abar`; non-periodic (outflow) derivatives (no circulant wrap,
so no inner box-mode artifact); validated line-Hilbert for `H(g)`.

## Spectral structure (analytic, then confirmed numerically)

In `z = ln(Y-1)` the transport `-2(Y-1)d_Y` is constant-coefficient `-2 d_z`, so the
operator has **continuous (essential) spectrum**, not just point eigenvalues. Its
real-part location is **norm-dependent** (the weight moves it; it is `<=0` and
pushable strictly left), so the stability verdict must rest on the **point**
spectrum. Sectors:
- **Omega-shape (`M_GG`)**: transport continuum (essential) + the nonlocal profile
  coupling. The question: does the coupling create a convergent, profile-localized
  point eigenvalue with `Re >= 0`?
- **Theta-shape (`M_hh`)**: pure transport (`c_l+2c_omega=0`, no damping) ->
  marginal essential spectrum at `Re=0` (the normalization / `nu` direction).
- **edge/`xi`**: appended symmetry mode at `lambda = c_l-1 = 1` (project out).

## Result (frozen mu=nu=0), validated

`validate/stability_edge_spectrum.py`:
- **`eig(M_GG)` with `Re>0`: 0** (with AND without the nonlocal term) — the nonlocal
  profile coupling adds no `Re>0` mode.
- **Converged, profile-localized point modes in `Re in (-3,2)`: NONE**, across
  `N:1600->2400`, `h0:2e-3->1e-3`, `L:160->320`. Every eigenvalue in the physical
  window DRIFTS under refinement => that part of the spectrum is purely
  essential/transport, not point spectrum.
- The **rightmost eigenvalues are far-field box modes** (peak at `x*=L`, low hi-k) —
  the truncated transport continuum at `Re<=0`, not profile modes.
- The predicted **edge/`xi` symmetry mode appears at `lambda = +1.000`** (validation
  target (i) met), and the **`Theta` sector is marginal at `Re=0`** as predicted.

**Verdict:** in the frozen sector the desingularized operator has **no genuine
unstable mode** — no point spectrum with `Re>0`, and the only `Re>=0` spectrum is
the two identified symmetry/normalization modes. This is **numerical evidence
consistent with Conjecture 2.4** (asymptotic stability of the singular profile),
obtained on the formulation that the smooth-profile CAP framework lacks.

## Honest status / not claimed

- **Frozen `mu=nu=0`.** The scaling/normalization modes (`mu=delta c_l`,
  `nu=delta c_omega`) are not yet added; they handle the marginal `Theta` direction
  and the scaling symmetries and must be projected out for the full asymptotic
  verdict. Structurally they append two more modes (like `xi`) with constraints.
- **Essential spectrum / spectral gap.** The verdict rests on the absence of point
  spectrum `Re>=0`; a publication-grade *spectral gap* needs a far-field weight (or
  mapped infinite domain) to push the transport continuum strictly left and remove
  the far-field box modes. The current near-edge nonlocal discretization also throws
  a deep-stable spurious cloud (`Re<-100`) that is harmless to the verdict but should
  be cleaned (cleaner near-edge quadrature, or the `z`-coordinate which makes the
  transport constant-coefficient and the essential spectrum a clean vertical line).
- **Not a proof.** Numerical, single-method, frozen. No nonlinear cross-check yet.
- **No interval arithmetic / Lean.** The CAP closure is a later, separate step.

## Nonlinear cross-check (gate (iii)) — edge mode CONFIRMED; shape test pending

`validate/nonlinear_crosscheck.py` (PASS):
- **DECISIVE.** The edge/`xi` mode `lambda = c_l-1 = 1` is confirmed by the EXACT
  nonlinear flow, grid-free: the edge-displaced singular profile `Omega_bar(X-s)`
  has closed-form velocity `U(X)=U_bar(X-s)-U_bar(-s)` (H commutes with translation),
  giving exact edge dynamics `ds/dtau = 2 sqrt(1+s)(sqrt(1+s)-1) -> s` as `s->0`
  (rate `1.0000`), and `X_e(tau)-1 = s0 e^{tau}` to machine precision. The naive
  no-feedback estimate would give rate `2`; the nonlinear flow gives `1`, confirming
  the `-xi` velocity self-feedback that is the key correction in the derivation.
- **LIMITATION.** Forward-marching a regularized singular profile (`c_l=2,
  c_omega=-1` fixed) on the equispaced grid DECAYS, and the decay depends on the
  regularization `delta` (ratios 0.67/0.64/0.52 for `delta=0.10/0.05/0.025`) -> a
  regularization/resolution artifact (the smoothed edge sits below the critical
  saddle amplitude), NOT a stability verdict. The shape-stability forward march
  needs the mapped-grid solver.

## Edge-RESOLVED confirmation in z = ln(Y-1) (cleaner, + a methodological catch)

Code: `src/stability_edge_z.py`, `validate/stability_edge_z_spectrum.py` (PASS).
Rebuilt the operator on a uniform `z=ln(Y-1)` grid (clusters at the edge, compresses
the `1/sqrt(Y)` tail). Transport `-abar G_Y = -2 d_z` is constant-coefficient. The
support-only Hilbert is done with a non-equispaced p.v. quadrature (subtraction
form), validated against `H((Y-1)^{-1/2})=0` on `Y>1` to ~1.5% (10x better than the
equispaced tool's ~15%; residual = tail truncation + resolution; verdict robust to it).

**Methodological catch (important).** The edge `z->-inf` is an INFLOW boundary
(characteristics move `+z`). Omitting the inflow/regularity BC, the discretization
re-admits the **inadmissible homogeneous edge branch** `G~(Y-1)^{-lambda/2}` (which
blows up at the edge for `Re lambda>0`) as a SPURIOUS point mode at `Re~+1.7` that
DRIFTS with the `z_min` cutoff (+1.70/+1.73/+1.75/+1.74 at `z_min=-12/-14/-16/-18`).
This is exactly the Robin/regularity edge BC flagged in the derivation. Imposing
`G(z_min)=0` removes it cleanly: `max Re eig(M_GG)` goes `+1.7 -> -0.001`, robust in
`z_min`. (A naive "report max Re" without the BC would have produced a
plausible-but-WRONG instability — the trap the handoff warned about.)

**Clean verdict (edge BC on, frozen mu=nu=0):** `eig(M_GG)` with `Re>0` = **0**;
rightmost `~ -0.001 -> 0^-` under refinement (the marginal transport continuum);
edge/`xi` mode appended at `lambda=1`; `Theta` marginal at `Re=0`. Confirms the
equispaced verdict, now **edge-resolved and clean** (no spurious near-edge cloud, no
box modes). The essential spectrum sits at `Re=0` in the bounded-`G` norm; a
far-field decay weight shifts it strictly left (asymptotic stability in that weighted
norm) — that is an analytic statement about the continuous operator (a finite-matrix
similarity is eigenvalue-preserving, so it does NOT show the gap numerically).

## Next (in priority order)

1. **Mapped-grid nonlinear shape test** — redo the forward march on the `z`-grid
   (edge now resolved) keeping the `1/sqrt(Y)` tail, to test shape stability
   nonlinearly (the piece the equispaced march can't reach).
2. **Add `mu,nu`** with the paper's normalization conditions; project out the
   symmetry subspace `{xi@1, mu, nu}`; read the asymptotic verdict on the remainder.
3. **`z=ln(Y-1)` coordinate + far-field weight** for a clean point/essential
   separation and a quantitative spectral gap (and to kill the spurious clouds).
