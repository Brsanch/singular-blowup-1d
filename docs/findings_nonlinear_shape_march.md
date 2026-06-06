# Finding: nonlinear shape march — the convective-weighted stability survives to
# FINITE amplitude (frozen base); no surprise nonlinear instability

Date 2026-06-06. Code: `validate/nonlinear_shape_march.py` (PASS), built on the
verified linear operator `src/stability_edge_z.build_blocks_z`. This is THE NEXT
CHIP #1: extend the linear convective spectral-gap result (`findings_spectral_gap.md`,
infinitesimal) to **finite amplitude**.

## The question

The frozen-sector LINEAR verdict is complete: no Re>0 point mode, convective gap
`2γ` in `L²(e^{-2γz}dz)`, symmetry modes weight-excluded, and no transient growth
for physical perturbations (`G=1`, `findings_transient_growth.md`). The linear
picture only governs *infinitesimal* perturbations. The genuine "no surprise
instabilities" claim is a statement about *finite-amplitude* dynamics: does a
finite perturbation of the singular profile still decay, or does the nonlinearity
surface an instability the spectrum cannot see?

## Method (reuse the verified operator; add only the one quadratic term)

Frozen base `Ω̄=1{Y>1}/√(Y-1)`, `Θ̄=(π/2)1{Y>1}`, `c_l=2, c_ω=-1`. Perturbation
`Ω=Ω̄+dΩ`, `dΩ=(Y-1)^{-1/2}G`, `Θ=Θ̄+h`. eq.(2.4) splits EXACTLY as
`dpert_τ = L(dpert) + N(dpert,dpert)`:

- `L` = the **verified** frozen linear operator (`M_GG, M_Gh, M_hh`) whose spectrum
  gives the gap. Reused unchanged — so the linear limit is an exact self-check.
- `N` = the only quadratic term, the perturbation self-advection `-dU·dpert_Y`
  (`dU = ∫₀^Y H(dΩ)` = the perturbation velocity). In the bounded G-variable:
  - `N_G = -(dU/(Y-1))·(G_z - G/2)`   [`= (Y-1)^{1/2}·(-dU·dΩ_Y)`]
  - `N_h = -(dU/(Y-1))·h_z`           [`= -dU·dΘ_Y`]
  with `dU/(Y-1) = Q_op@G`, the bounded edge form (exponential-averaging
  convolution, weights ≤1; the same conditioning fix as `findings_transient_growth.md`,
  added additively to `build_blocks_z`).

Non-moving frame (identical to `spectral_gap.py`); edge BC `G(zmin)=h(zmin)=0`;
bulk Gaussian bump IC (the `ξ`/symmetry mode is not excited and is weight-excluded).
RK4 time march; weighted decay rate `-d/dτ ln‖·‖_w` fit over the clean window.

## Validation ladder (passed before believing the verdict)

- **GATE A** (RK4 linear, N off) reproduces `spectral_gap.py` EXACTLY:
  rate `0.054 / 0.552 / 1.057 / 1.675` at `γ=0/0.25/0.5/0.75` (predicted `2γ`).
  → the time-stepper is correct.
- **GATE B** (N on, `A→0`) converges to the same linear rate (1.057 at γ=0.5 for
  `A≤0.01`). → the nonlinear term vanishes correctly as `O(A)`.

## Result — finite-amplitude sweep (γ=0.5 convective norm)

| A | rate(γ=0) | rate(γ=.25) | rate(γ=.5) | rate(γ=.75) | peak ‖·‖_w | edge hi-k |
|---|---|---|---|---|---|---|
| 1e-3 | +0.054 | +0.552 | +1.057 | +1.675 | 1.000 | 0.000 |
| 0.1  | +0.056 | +0.552 | +1.056 | +1.672 | 1.000 | 0.000 |
| 0.5  | +0.062 | +0.553 | +1.049 | +1.659 | 1.000 | 0.000 |
| 1.0  | +0.071 | +0.555 | +1.042 | +1.640 | 1.000 | 0.000 |
| 2.0  | +0.089 | +0.563 | +1.028 | +1.599 | 1.000 | 0.000 |
| 4.0  | +0.125 | +0.591 | +1.015 | +1.519 | 1.000 | 0.000 |

- The convective (`γ=0.5`) decay rate stays **strictly positive** up to `A=4`
  (4× the profile-scale perturbation), and `peak ‖·‖_w/‖·‖_0 = 1.000` — **monotone
  decay, no transient growth**, no grid-scale edge content (`hi-k = 0`).
- The nonlinearity is **active, not negligible**: it measurably *raises* the
  unweighted/low-`γ` rate (`0.054→0.125`: the self-advection adds extra outflow
  transport, sweeping energy to `+z` faster) while mildly *lowering* the high-`γ`
  weighted rate (`1.06→1.02`). Both directions stay decaying — the perturbation is
  carried out faster, not amplified.
- **Convergence** (`A=0.5, γ=0.5`): rate `+1.049` is IDENTICAL across
  `(N,zmin) = (1100,-12), (1400,-14), (1800,-16)`; peak `1.000`, hi-k `0`. The
  result is resolution-robust, not a discretization artifact.
- **Coupled (Ω,Θ)** (both bumps, `A=0.5`): a small bounded transient (peak `1.11`)
  then decays (`γ=0.5` rate `+0.935`). The Θ-coupling adds a modest non-normal
  bump, no instability.

## Verdict

The frozen-base **finite-amplitude** shape perturbation of the singular Hou–Luo
profile **decays in the convective weighted norm** for every tested amplitude
(`A` up to 4), with the nonlinear rate → the linear gap `2γ` as `A→0` and no
transient growth. The nonlinearity transports energy outward slightly faster; it
does not destabilize. **No surprise nonlinear instability.** This extends the linear
convective-gap verdict to the nonlinear, finite-amplitude regime — numerical
evidence consistent with Conjecture 2.4.

## Honest status / not claimed

- **Frozen base.** `Ω̄, Θ̄, c_l, c_ω` are held fixed; the perturbation self-advects
  but the base does not co-evolve and the exponents are not re-determined by a
  dynamic normalization. This tests "does a finite perturbation of the fixed profile
  decay", NOT "does the full nonlinear flow relax to the profile with self-adjusting
  `c_l, c_ω`." The latter (co-evolving base + dynamic normalization + moving edge) is
  the referee-facing gold standard and the recorded next obstruction.
- **Convective / weighted, not absolute.** Decay is in `L²(e^{-2γz}dz)`; the
  unweighted `L²(dY)` norm is marginal (the perturbation is swept to the far field
  without pointwise decay) — the standard sense for transport-dominated self-similar
  operators (Chen–Hou / De Huang lineage), exactly as in the linear result. NOT
  absolute decay.
- **Numerical, not a proof.** Finite edge-resolved discretization, ~1.5% Hilbert,
  RK4. No interval arithmetic / Lean. The grid-scale edge mode (the discretized
  `(Y-1)^{-1/2}` singularity) is unexcited here (hi-k=0) but is not function-
  analytically bounded.

## Next obstruction (handoff)

The remaining numerical step is the **co-evolving-base dynamic-rescaling march with
the paper's explicit `μ,ν` normalization** (next chip #2/#1-full): integrate eq.(2.4)
in the moving edge-fixed frame from a degenerate IC, letting `c_l, c_ω` adjust by the
normalization conditions, and confirm relaxation to `(Ω̄, c_l=2, c_ω=-1)` without
exponent drift. The known wall: the smooth/forward march is a saddle (dead-end #1)
and long-time marching fights the edge-resolution wall — it needs a stabilized
normalization (Newton/BVP-style or projected). After that, the precise weighted
space + interval-arithmetic/Lean (CAP) closure remains, and is scoop-risky.
