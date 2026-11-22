# Findings — hostile review of the global nonlinear weak-edge compatibility claims

**Date:** 2026-11-21 (review), added to the repo 2026-11-22
**Reviewed:** commit `12524e5` — `docs/review_packet_2026_11_21_weak_edge_nonlinear_manifolds.md`, the three
notes it lists (`findings_session_2026_11_21_nonlinear_cusp_manifold.md`,
`..._pure_step_nonlinear_manifold.md`, `..._finite_sobolev_pure_step.md`), the three gates (`NM-1..3`,
`PM-1..3`, `FS-1..3`), and `paper/stability_spectra_singular_houluo.tex` around `eq:dyn` and
`rem:weak_edge_jump`.
**Review form:** independent — the analytic audit was done by one Claude session that did not write those
notes (it is the author of LEM-10, which they build on); the numerical probes were run by a separate Claude
agent with fresh context, with its own real-space operator and its own solver copies. Not expert review.
The analytic findings were not cross-checked by a second reviewer.
**Scripts:** `validate/review_2026_11_21_weak_edge/` (the numerical agent's; kept unedited, not re-run by the
analytic reviewer; they reference a session scratch directory for outputs).
**Status of the reviewed material:** as of this note the three findings notes, their gates and the review
packet have NOT been changed in response. `rem:weak_edge_jump` in the manuscript has not been changed.
**Does not claim:** that weak-edge compatible solutions do not exist; that any theorem is false.

Severity: **S1** a reader could not tell it is wrong before relying on it · **S2** wrong or unsupported but
announces itself · **S3** wrong but harmless · **doc** wording.

## Verdict by claim (numbering of the review packet)

| Claim | Verdict |
|---|---|
| 1 — autonomous s-clock reduction | holds (re-derived by hand: the system on `v > 0`, dilation covariance of `Q`, the reduced system, the `r(s)` equation, the `V` formula; `V = r·a(s) + O(r³ log r)` near the horizon) |
| 2 — √v Θ-cusp strip-analytic manifold | not established: the contraction is false in the stated norm for the unmodified cusp (A1), and the note's own §5 lists its key estimate as open (A2) |
| 3 — pure-step strip-analytic manifold | the constructed solution survives independent real-space checks on a clean domain (D); "automatic compatibility" is by construction (A5, C1); the contraction is a sketch |
| 4 — finite-Sobolev `H³_{1/4}` closure | does not close as stated (A3); the phase-cancellation bound is false at `k = 0` (A4) |

## A. Genuine mathematical flaws / gaps

A1 [S2] Track A contraction is false in the stated norm for the unmodified cusp. With `h0 = 0` and forcing
`F(0,s) != 0`, `h = pi F(0)(1-e^{-us})/u` exactly for constant `F(0)`, so `W = s F(0) Psi(w/s^2)` and
`W^(m) = -2 cos(pi m) Gamma(2m+2) Gamma(2m-1) s^{1-2m} F(0)`. On any line `sigma2 > 1/2` the norm is
`~ s^{1-2 sigma2}` -> infinity as `s -> 0+` (s^{-1/4} at 5/8, s^{-1/2} at 3/4). The cone-norm hypothesis
`d ||W||_rho <= M_W` in `X_rho` fails at `s = 0+`. The headline "unmodified linear cusp, no tuning" and the
proof sketch need different things: either the modified cusp `h0 = pi c0(0)/(u+1)` (a G0-dependent O(eps^2)
correction of the initial cusp) or a time-weighted norm. Mild solution in `L^2_{1/4}` alone is fine.

A2 [S2] Track A is labelled "AFFIRMATIVE / Established" while its own §5 lists the `sigma2` pole-splitting
estimate for `S_g` and `T` as open. `S_g(t)` is NOT bounded on a line `sigma2 > 1/2` for data with
`a != 0`: `||S_g(t)U||_{sigma2} ~ |a(U)| t^{-(2 sigma2 - 1)}`; R->G class `~ t^{-2 sigma2}` (non-integrable).
The note's "0 derivative loss, uniformly in s" for `S_g` on `X_rho` is false on the second line. A weakly
singular Nishida scheme may close; none is written. Tier: sketch.

A3 [S2] Track C does not close as stated. `b = Q G~ = -2 N(u h)` gives `b in H^j(R)` (unweighted) iff
`h in H^{j+3/2}_{1/2}`. The energy gives `h in L^inf_s H^3_{1/4}` and the dissipation `L^2_s H^3_{1/2}`.
The symmetric-part bound needs at least `int |eta||b^|` i.e. `b in H^{3/2+delta}` i.e. `h in H^{3+delta}_{1/2}`:
short by `delta > 0` derivatives at the edge side. `h in H^3_{1/4}` controls only `e^{z/4} b in H^2`.
The parenthetical "(z<=0 bound supplied by u-dissipation / H^4_{3/4} after splitting the rank-1 residue)" is
not an argument; on a line > 1/2 the conjugation is not exact (rank-one defect). Claim 4 is unproven.

A4 [S1] FS-1 bound is false at k = 0 (and unproven at k = 1). Kernel bracket with the true amplitude behaves
like `|eta|^{1+q}`, `q = 2-k`, when the input frequency is near 0. Numerics: for `b` barely in `H^2`,
`||sym_0||` = 12.9 / 24.7 / 47.4 / 91.5 / 176.6 at cutoff 16..257; narrow Gaussians: ratio to `||b||_{H^2}`
doubles per halving of width. k = 2, 3 and the single order-3 weight hold. Harmless for the k = 3 argument,
but the note says "for k in {0,1,2,3}" and the gate prints "H^0: 3.6546" as verified.

A5 [S2] Selection / non-uniqueness is never stated. On `v > 0` the linearized operator annihilates constants
(`Q(1) = 0`), so for ANY edge-amplitude function `alpha(tau)` there is a solution from zero data
(`G = alpha + G1`, `G1 = -int S(tau-tau') alpha'(tau') 1 dtau'`). The edge is an inflow boundary; the edge ODE
is satisfied by all of these. `G(0,tau) = 0` is a boundary condition imposed by working in `L^2_{1/4}`;
"automatic compatibility" is by construction. Numerically confirmed: `a = 0` and `c1 = 3 a(forcing)` hold
for wrong or unrelated nonlinearities. The solution regular up to the edge from the same data violates
compatibility; uniqueness of weak solutions is not addressed.

A6 [S2, manuscript] `rem:weak_edge_jump` is now contradicted: "H[g](1) = 0 forces G = 0 at linear order" and
"no nontrivial solution of the linearized flow preserves B = A^2" hold only in the edge-analytic class.
Counterexample: `S(tau) N k0`, `k0 in C_c^inf` (LEM-10).

## B. Repairable proof gaps

B1 `||QG||_inf` justification ("Cauchy–Schwarz on Re zeta = ±1/4") is wrong when `a != 0` (pole of cot at 0;
`QG` is not in `L^2_{1/4}`); the bound is true via `QG = a + L^2_{1/4}` part.
B2 Lipschitz of `c0(s)`: looks fine — `a(P(G))` takes no derivative after one integration by parts and needs
no analyticity radius — but is not written.
B3 Weak-edge equivalence across X = 1: no hidden delta or (X-1)^{-1/2} defect found given `G(0,tau) = 0`,
`eta(0+,tau) = 0`; conditional on the paper's convention that `U_X Omega` carries `-(pi/2)A^2 delta`
(boundary value of the analytic square; the pointwise a.e. product is 0). That convention is unchecked
against arXiv:2604.01868.
B4 Horizon: `V = r a(s) + O(r^3 log r)`, Lipschitz, transversal crossing — fine; but no gate tests it (C3).

## C. Numerical / gate overclaims

C1 [S1] Every "a = 0" metric (PM-2, PM-3, FS-3, NM-2 edge balance) is true by construction; PM-2's `max|a|`
fails its own 5e-6 threshold under dt refinement (2.8e-6 -> 1.5e-5).
C2 [S1] PM-2 and NM-3 printed residuals are periodic-wraparound values at the mask edge (6.7e-5 -> 8.9e-6 on
[-48,48]; NM-3 9.9e-6 -> 4.6e-8 on [-64,64]). In the gate's PM-2 solution max|G| on z in (12,38) is 0.30.
C3 [S1] PM-3 / NM-3 physical-time checks pass with V_G replaced by -V_G, 3 V_G or 0.05 sin(7r): they cannot
see V_G. The V_G formula, the gauge, s_inf and tau(s) are untested by any gate.
C4 [S1] FS-3 H^3 diverges with Nz on the gate domain (2.7e-2, 5.6e-2, 0.29, 1.6); converged 2.5767e-2 on
[-72,72]. Attribution to wraparound confirmed; no interior roughness growth. The gate prints only L2/H1.
C5 [S2] The scheme goes NaN on long-left domains with the filter on; unfiltered and fine-grid runs are stable.
No masked PDE instability found.
C6 [S1, repo] `src/stability_edge_z.py` `Q_op` has no tail below z_min: its edge value is ~0 for any G.
C7 [doc] NM-1 "edge err 0.0" and NM-2(c) are tautologies.

## D. Survives

- Claim 1 algebra (hand re-derivation of (*), dilation covariance, (dagger), r equation, V formula).
- Equivalence (ddagger); `f -> pi F(0)`; `N'(1/2) = -2`; `c1 = 3 a(P)`.
- Phase-cancellation identity; top-order (k = 2, 3) symmetric-part bound; FS-2 (exact infimum 0.893).
- Conjugation identity on evolving profiles: real-space Q vs conjugated Q 5e-9.
- On clean domains the conjugated pure-step solution satisfies the real-space PDE to 3e-9 in the core
  (~1e-6 of the nonlinear term); `a` by independent tail-corrected quadrature ~1e-12; `c1 = 3a(P)` to 6 digits.
- Physical-time residual by an independent path: 4.3e-6 / 2.6e-6 at tau = 0.2 / 0.6 (0.6% / 2.5% of the
  nonlinear terms); tau = 1.5, 3.0 uninformative.

## Not checked

Independent real-space time integration; the coupled Track A system in real space; any written Nishida
contraction; Lean; the source paper's weak-solution convention.

## E. Numerical detail (from the independent agent)

**`a ≡ 0` holds for any forcing** (gate metric `max|Q_conj(G)(z=−15)|`, threshold 5e-6; `a_quad` = direct
quadrature with the `c₁w` tail; clean grid `[−80,80]`, no filter):

| Nonlinearity in the conjugated solver | gate metric on the gate grid | `a_quad`, clean grid | `G/w` at z = −14 vs `3a(nl)` |
|---|---|---|---|
| true `𝒫` | 1.8e-07 | −1.3e-13 | −1.4092e-03 vs −1.4091e-03 |
| advection sign flipped | 3.7e-07 | +3.0e-13 | +2.7657e-03 vs +2.7657e-03 |
| `G²` (no `Q` at all) | 2.6e-08 | −2.8e-14 | −3.2573e-04 vs −3.2577e-04 |
| `0.3·exp(−(z+3)²/2)·G` | 2.2e-07 | −4.2e-13 | −8.8426e-03 vs −8.8435e-03 |

**PM-2's `max|a|` is not converged in the time step** (clean domain): 2.79e-06 / 1.31e-05 / 1.45e-05 / 1.49e-05
at `N_s` = 240 (gate) / 480 / 960 / 1920, against its own threshold 5e-6.

**Printed residuals are wraparound:** PM-2 6.66e-05 (at s = 0.05, z = 11.97, the mask edge) → 8.85e-06 on
`[−48,48]`; NM-3 `res_G` 9.93e-06 → 8.75e-08 on `[−48,48]`, 4.63e-08 on `[−64,64]`. In the gate's PM-2
solution `max|G|` over z ∈ (12, 38) is 0.30 (15× ε).

**Physical-time checks cannot see `V_G`:** PM-3 and NM-3c PASS with `V_G` replaced by `3·V_G`, `−V_G` or
`0.05·sin(7r)` (residuals 3.3e-06 … 4.0e-06); they fail only through the `|s_∞ − 1|` bound when `V_G = 0.3·r`.

**FS-3 norms of `v = e^{−z/4}h` at s = 1.2** (gate domain `[−36,36]`, filter 2/3):

| `N_z` | L² | H¹ | H² | H³ |
|---|---|---|---|---|
| 1024 | 1.70278e-03 | 2.93548e-03 | 6.991e-03 | 2.716e-02 |
| 2048 | 1.70261e-03 | 2.93520e-03 | 7.020e-03 | 5.629e-02 |
| 4096 | 1.70267e-03 | 2.93532e-03 | 7.257e-03 | 2.863e-01 |
| 8192 | 1.70267e-03 | 2.93537e-03 | 8.905e-03 | 1.603 |

On `[−72,72]` and `[−96,96]` (matched `dz`): H² = 6.98623e-03, H³ = 2.5767e-02, monotone decreasing in s,
filter-independent. The H³ drift is periodic wraparound; no interior roughness growth was found.

**FS-1 symmetric-part norm for `b` barely in `H²`:** k = 0: 12.9 / 24.7 / 47.4 / 91.5 / 176.6 at cutoff
16 … 257 (unbounded); k = 1: 2.12 → 2.84 (slowly growing); k = 2: 0.790 (flat); k = 3: 0.498 (flat).
FS-2's exact infimum over all `v` is 0.8935 (> 0.5).

**Independent real-space residual of the pure-step equation** (own p.v. `Q`, exact `G_s`, clean
`[−104,24]`, no filter): 3.1e-12 (edge), 2.7e-09 (core), 2.0e-08 (bulk), 2.2e-07 (far); dropping the
nonlinear terms gives 2.2e-03 in the core. Physical-time equation by an independent path: 4.3e-06 / 2.6e-06
at τ = 0.2 / 0.6; τ = 1.5 and 3.0 are uninformative (finite-difference limited).

**Solver stability:** the filtered ETD-RK2 and Picard schemes go NaN on long-left domains
(`[−104,24]`; `[−72,72]` at `N_z` ≤ 4096), amplitude-dependent; unfiltered and fine-grid runs are stable and
agree to 6 digits.

**`src/stability_edge_z.py` `Q_op`** had no tail below `z_min` (edge value ≈ 0 for any `G`). That finding led
to `docs/findings_session_2026_11_22_qop_edge_tail.md`; the default operator has since been corrected.

## F. What a repair would need (reviewer's reading)

1. Track A: either adopt the modified cusp `h₀ = πc₀(0)/(u+1)` and say so in the statement, or use a
   time-weighted norm on the second line; write the weakly singular Nishida argument for `S_g` and `𝒯`.
2. Track C: control of `b` in `H^{3/2+δ}` at the edge side — a second weighted line with the rank-one
   defect handled explicitly, or a sharper commutator bound; restate FS-1 for `k ∈ {2, 3}` only.
3. State the inflow-boundary non-uniqueness (A5) and that `G(0,τ) = 0` is imposed.
4. Gates: domain of at least `[−64,64]` (or `[−104,24]` unfiltered) with per-region residuals and asserted
   `N_z` / `N_s` / domain convergence; tail-corrected quadrature for `a` with a mutation control; a
   physical-time check that depends on `V_G`; H² and H³ in FS-3.
5. Correct `rem:weak_edge_jump` (A6) and check the singular-product convention against the source paper (B3).
