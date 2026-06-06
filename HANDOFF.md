# HANDOFF — singular-blowup-1d

Read this + `README.md` + `docs/findings_*.md` first. The next step is named at
the bottom ("THE NEXT CHIP"). Repo: github.com/Brsanch/singular-blowup-1d (private).
Started 2026-06-05. Single-researcher + AI; tight measure→mechanism→remeasure loop.

## The bet (one paragraph)

Certify *singular* self-similar finite-time blowup of 1D nonlocal active-scalar
models (CLM / De Gregorio / gCLM / **Hou–Luo**). This INVERTS the NS/SQG
regularity wall: instead of ruling out blowup on the true equation (→ Type-I /
KNSS-Liouville rigidity, equivalent-to-conclusion), we *certify* blowup on toy
models where blowup is the answer. 1D ⇒ extreme resolution cheap; the profile
ODEs are the cleanest CAP/Lean targets. Chosen over SQG-patch / OSW / gCLM-
classification by a 27-candidate adversarial scout (`docs/candidate-scan.md`).

## Verified landscape (the two anchor papers, both real, checked)

- **arXiv:2603.25104** (Huang, Tong, Wang, Mar 2026): gCLM family `ω_t+a u ω_x=u_x ω`,
  singular profiles for `a<0` (explicit), `a=0` (CLM) rigorous, `a>0` regular.
- **arXiv:2604.01868** (Chen, **De Huang**, Li, Apr 2026): 1D Hou–Luo + 2D
  Boussinesq, two-stage `L∞→Lᵖ` blowup, **unbounded `(X−1)^{−1/2}` profile**,
  exponents `c_l=2, c_ω=−1, γ=2`. NUMERICAL ONLY — no rigorous existence, **no
  published stability spectra**. ← our hard-to-scoop target.

**Scoop risk:** De Huang's group owns the CAP machinery and is driving this hard.
Robust contributions: independent confirmation (done) + stability spectra (open).

**Target (Bryan picked):** "confirm + stability spectra, Hou–Luo lead."

## The Hou–Luo model + the certified base state

Model (`x∈ℝ`): `ω_t + u ω_x = θ_x`, `θ_t + u θ_x = 0`, `u_x = H[ω]`.
Dynamic-rescaling system eq.(2.4) (steady state = the profile; linearization =
the stability operator):
```
Ω_τ + (U + c_l X)Ω_X = c_ω Ω + Θ_X
Θ_τ + (U + c_l X)Θ_X = (c_l + 2c_ω)Θ,   U_X = H(Ω), U(0)=0
```
**Certified singular steady state (Thm 2.3), closed form:**
```
Ω̄(X) = 1{X>1}/√(X−1),  Θ̄(X) = (π/2) 1{X>1},  c_l=2, c_ω=−1,
Ū(X)  = 2√(1−X)−2  (X<1),  −2  (X>1)   ⇒  U(1)=−2 exactly,
ā(X) := Ū+2X = 2(X−1) on X>1.
```

## State of play (all committed + pushed)

| result | status | where |
|---|---|---|
| gCLM engine vs CLM exact | ✅ 2e-10 | `validate/validate_clm.py` PASS |
| Hou–Luo physical DNS `c_ω=−1` | ✅ <0.4% | `findings_houluo_physical.md` |
| Thm 2.3 singular profile certified (`U(1)=−2`, residual 2e-16) | ✅ | `validate/verify_singular_steady.py` PASS, `findings_singular_steady_verified.md` |
| dynamic-rescaling solver + "profiles are saddles" | ✅ | `findings_dynamic_rescaling.md` |
| stability spectrum (frozen `μ=ν=0`): resolved sector stable; leading mode = free-boundary artifact | ✅ evidence-for-stability | `findings_stability_spectrum.md` |
| edge-condition derivation (the analytic core) — done + cross-checked 4 ways + verified | ✅ | `findings_edge_condition_derivation.md`, `validate/probe_edge_condition.py` PASS |
| **desingularized free-boundary operator BUILT; frozen spectrum has NO genuine unstable mode** | ✅ evidence-for-Conj-2.4 | `src/stability_edge.py`, `validate/stability_edge_spectrum.py` PASS, `findings_edge_stability_operator_built.md` |
| nonlinear cross-check: edge mode `λ=c_l−1=1` confirmed by EXACT flow (rate 1 not naive 2); shape-march inconclusive (regularization wall) | ✅ edge-mode / ⚠️ shape pending | `validate/nonlinear_crosscheck.py` PASS |
| edge-RESOLVED operator in `z=ln(Y−1)`: clean spectrum, NO Re>0; caught+fixed spurious edge mode (inflow BC `G(z_min)=0` essential) | ✅ confirms verdict, edge-resolved | `src/stability_edge_z.py`, `validate/stability_edge_z_spectrum.py` PASS |
| **SPECTRAL GAP exists**: Ω-shape asymptotically stable in convective weighted norm `L²(e^{−2γz}dz)`, gap `2γ` (decay rate matches to ~5%); unweighted marginal | ✅ resolves asymptotic-vs-marginal crux (frozen Ω-sector) | `validate/spectral_gap.py` PASS, `findings_spectral_gap.md` |
| **Θ + μ,ν folded in**: weighted space IS the symmetry quotient (ν=constant & dilation/edge=ξ@1 excluded by weight); full coupled (G,h) gapped, no Re>0 | ✅ complete frozen-sector verdict | `validate/theta_munu_gap.py` PASS, `findings_theta_munu_gap.md` |
| **operator = correct linearization of trusted eq.(2.4) RHS** (rel.err ∝ ε→0); caught+fixed an upwind/central base-deriv inconsistency; weighted decay tendency <0 | ✅ nonlinear cross-check (decisive form) | `validate/nonlinear_linearization_check.py` PASS, `findings_nonlinear_linearization_check.md` |
| transient growth / non-normality SETTLED: smooth physical perturbations have G=1 (NO transient growth) at all resolutions; worst-case creep is a grid-scale edge mode (hi-k=1.0), not physical | ✅ linear stability robust; nonlinear supported | `validate/transient_growth.py` PASS, `findings_transient_growth.md` |
| **NONLINEAR shape march (finite amplitude, frozen base)**: convective weighted decay survives to `A=4` (4× profile scale); rate→linear gap `2γ` as `A→0` (GATE B), peak `‖·‖_w=1.000` (monotone, no transient growth), resolution-robust, edge hi-k=0; coupled (Ω,Θ) adds only a bounded `1.11` transient | ✅ extends gap to nonlinear; no surprise instability | `validate/nonlinear_shape_march.py` PASS, `findings_nonlinear_shape_march.md` |
| **co-evolving march STARTED**: dynamic normalization (`c_l=-U(1)` edge-stagnation + `c_ω` amplitude-pin) VALIDATED — profile is an exact fixed point (`c_l=2.0000000`, `c_ω=-1.0000000`, `dU(1)` functional rel.err 7.6e-8); fixed-Y dilation gauge term quantified to diverge as `e^{-z}` at edge ⇒ MOVING frame mandatory | ✅ machinery built; ⚠️ wall = needs moving frame (next) | `validate/coevolve_march.py` PASS, `findings_coevolve_march_started.md` |
| **co-evolving march MOVING frame**: λ=1 pure-ξ symmetry mode confirmed in the coupled (g,h,ξ) operator (rest Re≤0); along the stable trajectory the diagnostic exponents converge `c_l→2, c_ω→-1` with NO drift (resolution-robust), edge pinned, `‖g‖_w` decays at the gap rate; explicit forward forcing of the exponents is saddle-unstable (= dead-end #1) ⇒ needs PROJECTED normalization | ✅ λ=1 gate + diagnostic convergence; ⚠️ self-consistent forcing needs projection | `validate/coevolve_moving_frame.py` PASS, `findings_coevolve_moving_frame.md` |
| **co-evolving march SELF-CONSISTENT (relaxes!)**: runaway localized to the Θ-amplitude h-forcing; Θ-balance normalization `c_ω=-c_l/2` removes it (`c_l+2c_ω≡0`); FULLY self-consistent march (both `c_l,c_ω` fed back) RELAXES — perturbed profile → `(Ω̄,c_l=2,c_ω=-1)`, exponents converge to 6 digits, NO drift, resolution-robust, edge pinned ⇒ singular profile is a NONLINEAR ATTRACTOR of the dynamic-rescaling flow | ✅ co-evolving relaxation (removes frozen-base + diagnostic-only caveats) | `validate/coevolve_projected.py` PASS, `findings_coevolve_self_consistent.md` |

## File map

```
src/gclm.py            gCLM/De Gregorio engine (CLM-exact-validated)
src/clm_exact.py       CLM closed form (self-contained Riccati derivation)
src/houluo.py          physical-space Hou-Luo solver (c_omega=-1)
src/rescaled.py        R-line Hilbert (validated 1.4e-6) + singular profile + steady residual
src/dynamic_rescaling.py  eq (2.4) time solver (saddle finding)
src/stability.py       2N x 2N linearized operator (frozen mu=nu=0), closed-form U_bar
validate/              all gates + probes + the spectrum convergence study
docs/                  candidate-scan + 4 findings docs
```
Run any gate: `python3 validate/<name>.py`. numpy/scipy present. Pure-Python, 1D,
cheap (no GPU). No kernel-panic risk (no du/find/lake here).

## Dead ends (do NOT redo)

1. Forward time-marching of eq (2.4) to reach SMOOTH profiles — they are saddles
   (nonlinear via `U=H(Ω)`); constraint-`c_ω` runs away (+39), fixed `c_ω=−1`
   decays to 0. Use Newton/BVP if smooth profiles are needed (they are NOT needed
   for the singular-profile stability deliverable — singular profile is closed-form).
2. NAIVE regularized eigenvalues of the singular operator — leading modes are
   δ-divergent grid-scale EDGE artifacts (`−0.13 → +24+36i` as δ:0.14→0.035,
   `x*≈1`, hi-k~0.45). Do NOT report these as instabilities.
3. `c_l` from physical DNS — the self-similar core `(T−t)^{c_l}`, `c_l>1`, goes
   sub-grid immediately. Use dynamic-rescaling coords.

## DONE (2026-06-05): edge derivation + operator built + frozen verdict

Read `docs/findings_edge_condition_derivation.md` then
`docs/findings_edge_stability_operator_built.md`. The operator is **built**
(`src/stability_edge.py`) and validated (`validate/stability_edge_spectrum.py`):
frozen `μ=ν=0` spectrum has **NO genuine unstable mode** (0 eigenvalues with Re>0;
no convergent profile-localized point mode; rightmost are far-field box modes =
transport continuum ≤0; edge/`ξ` symmetry mode at `λ=1` as predicted; `Θ` marginal
at Re=0). **Numerical evidence consistent with Conjecture 2.4 in the frozen sector.**
NOT a proof; frozen; no nonlinear cross-check yet.

## THE NEXT CHIP — complete the stability verdict

FROZEN-SECTOR VERDICT COMPLETE, NOW INCLUDING FINITE-AMPLITUDE NONLINEAR: full
coupled (Ω,Θ) shape operator, edge-resolved, no Re>0 point mode, spectral GAP `2γ`
in the convective weighted norm; scaling symmetry modes (ν=constant, dilation/edge=
ξ@1) are exactly the weight-excluded directions; AND the finite-amplitude frozen-base
shape perturbation decays in that norm to `A=4` with no transient growth (`nonlinear_
shape_march.py`). Numerical evidence consistent with Conjecture 2.4. Remaining,
priority order:
1. ✅ DONE (2026-06-06): nonlinear shape march on the mapped grid, FROZEN base —
   finite-amplitude convective-weighted decay confirmed, no surprise instability
   (`findings_nonlinear_shape_march.md`). The frozen-base caveat is the open part →
   item 2.
2. **CO-EVOLVING-base dynamic-rescaling march with explicit `μ,ν` normalization**
   (the referee-facing gold standard; subsumes the old "explicit μ,ν" item). STARTED
   2026-06-06 (`findings_coevolve_march_started.md`): the dynamic normalization is
   built + VALIDATED (`c_l=-U(1)` edge-stagnation + `c_ω` amplitude-pin → profile is an
   exact fixed point, `c_l=2.0000000/c_ω=-1.0000000`; the `dU(1)=∫_0^1 H(δΩ)` functional
   matches direct quadrature to 7.6e-8). WALL QUANTIFIED: in the fixed-Y frame the `c_l`
   (dilation) gauge term `-δc_l·X·Ω̄_X` diverges as `e^{-z}` at the edge (482→3563→26330→
   194500 as zmin -10→-16) because `X·Ω̄_X ~ (Y-1)^{-3/2}` is more singular than the
   profile. ⇒ the MOVING edge-fixed frame is MANDATORY (the frame `M_GG` already uses).
   ✅ DONE 2026-06-06 (`findings_coevolve_self_consistent.md`, `validate/coevolve_projected.py`):
   the FULLY self-consistent co-evolving march RELAXES. The runaway was localized to the
   Θ-amplitude h-forcing (not the Ω-amplitude); the Θ-balance normalization `c_ω=-c_l/2`
   makes `c_l+2c_ω≡0` so that forcing vanishes. With both `c_l=-U(1)` and `c_ω=-c_l/2`
   co-evolving and the edge pinned (`μ=-V_g(1)`), a perturbed profile relaxes to
   `(Ω̄,c_l=2,c_ω=-1)`: exponents converge to 6 digits, NO drift, resolution-robust across
   3 grids ⇒ the singular profile is a NONLINEAR ATTRACTOR of the dynamic-rescaling flow.
   (Diagnostic + moving-frame λ=1 gate also stand.) Coupled (Ω,Θ)-shape perturbation also
   relaxes (no coupling-driven instability); stress test (amplitude≤2, shapes: bulk/near-
   edge/far/high-k/two-bump/Θ-dominant) finds NO counterexample. ✅ CHIP #2 COMPLETE
   2026-06-06: paper-match done (`findings_paper_normalization_match.md`) — both anchor
   papers (2604.01868, 2603.25104) verified real; our `c_l=-U(1)` IS the paper's eq (3.2)
   condition 1 (holds identically), the profile satisfies condition 2 (`U_x(0)=-1`) and
   exponents (2,-1), and along our march `U_x(0)→-1` (the paper's `c_omega` condition,
   satisfied asymptotically; our Θ-balance `c_omega` agrees at/near the profile).
3. **Pin the precise weighted space / max gap; rigor.** Then interval-arithmetic/Lean
   closure (the CAP step; scoop-risky, hard). NB transient growth (physical G=1) and
   finite-amplitude nonlinear decay are now settled; rigor + co-evolving base are the
   remaining frontier.

⚠️ LESSON (this session): on the `z`-grid the edge is an INFLOW boundary; WITHOUT the
inflow/regularity BC `G(z_min)=0` the discretization shows a SPURIOUS unstable edge
mode at `Re~+1.7` that drifts with `z_min`. Always impose it; always check `z_min`/N
drift before believing any Re>0 mode.

---
### Reference: corrected formulation (built in stability_edge.py / stability_edge_z.py)

- **Moving (edge-fixed) frame is mandatory.** `Y = X − εξ(τ)`, keep `Ū(Y−εξ)`
  EXACT. The fixed-X-frame `−ξΩ̄'` gives a SPURIOUS infinite edge velocity
  (`xi(1-Y)^{-1/2}`, Taylor breakdown at the vertical-tangent edge) — that IS the
  recorded "δ-divergent edge artifact." In the moving frame the perturbation
  velocity is finite: `V_Y=H(g)`, `V(0)=−ξ` ⇒ `V(Y)=V_g(Y)−ξ`, `V_g(Y)=∫_0^Y H(g)`.
- **Edge ODE (derived 4 independent ways, all agree):**
  `ξ̇ = (c_l−1)ξ + V_g(1) + μ`,  `V_g(1)=∫_0^1 H(g)dY`. With `c_l=2`: `ξ̇=ξ+V_g(1)+μ`.
  The coupling is to the VELOCITY `δU(1⁺)=V_g(1)−ξ`, NOT `Hω(1)`(=gradient).
- **Shape eqs (Y>1, frozen μ=ν=0 core):**
  `g_τ = −ā g_Y + c_ω g + h_Y − [∫_1^Y H(g)] Ω̄_Y`,  `h_τ = −ā h_Y`. The nonlocal
  term is regular (`~ ½H(g)(1)(Y−1)^{−1/2}` at edge). `h` is a slaved transport mode.
- **Edge BC = Robin, not Dirichlet.** Indicial: `g~(Y−1)^{−(λ+1)/2}`; weighted
  space admits the homogeneous branch only for `Re λ≤0`. Any `Re λ>0` mode must
  DROP the homogeneous branch and keep only the forced `(Y−1)^{−1/2}` part with
  edge coeff `B=H(g)(1)/(2λ)`. Impose THAT, not `Hω(1)=0`.
- **Use closed-form `Ū, Ū_X, Ū_XX`** near the edge; never discretely Hilbert the
  singular `Ω̄'`. **Far field: NON-periodic (outflow) derivative or a mapped
  infinite domain — NOT a circulant `D`** (the circulant periodic-wrap fighting the
  line-Hilbert is the source of the spurious box modes at x*=40; commutator O(1)
  there, machine-zero at the edge — see probe (B)).
- **Symmetry modes to project out:** the pure-`ξ` mode is **predicted at
  `λ = c_l−1 = 1`** (time-translation-type); plus the `μ,ν` (scaling/amplitude)
  directions. Only after projecting these is the sign of `max Re λ` the verdict.

- **Validation gates:** (i) the `λ=1` pure-`ξ` symmetry mode appears (concrete
  falsifiable check); (ii) eigenvalues converge under mesh + map-parameter
  refinement, are edge-resolved (low hi-k), with NO box modes; (iii) cross-check
  the verdict against a NONLINEAR dynamic-rescaling DNS (`src/dynamic_rescaling.py`)
  from a DEGENERATE IC (`ω_x(0)=θ_xx(0)=0`) with a stabilized normalization. HARD
  STOP if exponents drift under refinement.

⚠️ This is the one place a careless discretization yields a plausible-but-WRONG
spectrum. Validate every step; do not report an unstable λ whose eigenmode is
edge-spiked/grid-scale, δ-dependent, or a box mode.

## Bryan's working preferences for this repo
Single-agent serial; local-verify-primary; commit each gated milestone; never
report a spectrum/instability without a convergence + resolution check
(no-overselling); state the obstruction at maximum resolution and hand off the
non-interpolative analytical leap rather than grinding past it.
