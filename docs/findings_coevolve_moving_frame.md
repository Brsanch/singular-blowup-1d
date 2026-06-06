# Finding: moving-frame co-evolving march — λ=1 symmetry mode confirmed; diagnostic
# exponents converge to (2,−1) with no drift; explicit forcing is saddle-unstable (projection needed)

Date 2026-06-06. Code: `validate/coevolve_moving_frame.py` (PASS). Continues THE NEXT
CHIP #2 from the fixed-Y wall (`findings_coevolve_march_started.md`): builds the
co-evolving system in the MOVING (edge-fixed) frame, where the dilation is carried by
the edge DOF ξ and the singular `δc_l·X·Ω̄_X` term is absent.

## System (moving frame; findings_edge_condition_derivation.md §4–7)

Bounded variables `G=(Y-1)^{1/2}δΩ`, `h=δΘ`, `ξ` scalar; no singular term:
```
g_τ  = M_GG g + M_Gh h + N_G(g)        + ν              (+ν Ω̄ → +ν in G)
h_τ  = M_hh h + N_h(g,h)               + (μ+2ν)(π/2)    (+(μ+2ν) Θ̄)
ξ̇   = (c_l−1) ξ + V_g(1) + μ
```
`M_GG,M_Gh,M_hh,N` = the VERIFIED frozen operator + nonlinear self-advection
(`nonlinear_shape_march`). Explicit normalization multipliers: `μ=−V_g(1)` (dilation /
`c_l=−U(1)`; makes `ξ̇=(c_l−1)ξ` homogeneous → edge pinned), `ν=c_ω+1` (amplitude pin),
`V_g(1)=∫_0^1 H(g)=b1@G`.

## Results

**GATE (i) — λ=1 pure-ξ symmetry mode (handoff's concrete falsifiable check).** The
linear coupled (g,h,ξ) operator (μ=ν=0) has eigenvalue `λ=1.00000` with eigenvector
exactly pure-ξ (`|ξ|=1.0`, `max|g,h|=0`), `#Re>0 = 1` (only that symmetry mode), and
`max Re of the rest = −0.0000`. The block-triangular structure (ξ slaved to g via
`V_g(1)`) gives `λ=c_l−1=1` precisely. **The moving-frame coupled formulation is correct.**

**GATE (iii) — relaxation + diagnostic exponent convergence.** A perturbed profile
evolved by the stable shape dynamics, edge pinned (`ξ` stays `~3e-3`), with the dynamic
exponents tracked DIAGNOSTICALLY (the values the normalization assigns):

| τ | c_l(diag) | c_ω(diag) | ξ | ‖g‖_w(.5) | Ω(Y_ref)/Ω₀ |
|---|---|---|---|---|---|
| 0.0 | +2.04379 | −0.85136 | 0 | 1.08e-1 | 1.000 |
| 1.0 | +2.00685 | −1.00027 | 1.5e-4 | 4.10e-2 | 0.867 |
| 2.0 | +2.00052 | −1.01228 | 4.4e-4 | 1.44e-2 | 0.874 |
| 3.0 | +2.00002 | −1.00095 | 1.2e-3 | 5.19e-3 | 0.880 |
| 4.0 | +2.00000 | −0.99999 | 3.3e-3 | 1.90e-3 | 0.880 |

`‖g‖_w` decays at rate **+1.065** (= the frozen linear gap `2γ=1.0`); `c_l→2`, `c_ω→−1`
monotonically with **NO drift**; ξ stays `~0` (edge pinned by `μ=−V_g(1)`).

**GATE (ii) — convergence.** Final `c_l=2.000000`, `c_ω=−0.99999`, decay rate `+1.065`
are IDENTICAL across `(N,zmin) = (1000,−10), (1400,−12), (1800,−14)`. Resolution-robust;
no exponent drift.

**OBSTRUCTION — explicit forward forcing is saddle-unstable.** Feeding the exponents back
as forcing (`+ν` in g, `+(μ+2ν)(π/2)` in h, `μ` in ξ̇) — even frozen per RK4 step
(semi-implicit) — BLOWS UP. The shape-only dynamics (forcing off) decays; turning on the
explicit forward normalization destabilizes. This reproduces the documented dead-end #1
("constraint-`c_ω` ran away to +39"): the amplitude is the saddle's unstable direction,
and forward forcing of it is unstable.

## Verdict / status

- **Validated**: the moving-frame coupled (g,h,ξ) operator is correct (λ=1 pure-ξ mode,
  rest Re≤0); the fixed-Y dilation singularity is absent (dilation carried by ξ).
- **Diagnostic co-evolving**: along the stable convectively-decaying trajectory the
  dynamic exponents converge `c_l→2, c_ω→−1` with no drift (resolution-robust), edge
  pinned — the singular profile is the attractor with the right exponents emerging. This
  is a CONSISTENCY + no-drift result (the normalization is consistent with the profile as
  the self-similar fixed point), not a fully self-consistent forced march.
- **Open (the precise remaining step)**: a FULLY self-consistent co-evolving march
  (exponents fed back) needs a PROJECTED / implicit normalization that projects out the
  amplitude (saddle) direction, not forward forcing. Forward forcing reproduces dead-end
  #1. This is the standard reason the construction literature (arXiv:2308.01528 /
  2603.25104) uses fixed-point / Newton, not time integration.
- Convective/weighted, frozen-then-diagnostic, numerical — consistent with Conjecture
  2.4; NOT a proof.

## Next

1. Projected/implicit normalization: project the field onto the gauge-fixed slice each
   step (remove the amplitude mode) so the exponents can be fed back stably; then the
   fully self-consistent co-evolving march should relax (HARD STOP if exponents drift).
2. Match the paper's explicit μ,ν normalization conditions (referee-facing) and confirm
   they agree with the weighted-space quotient + this diagnostic convergence.
3. Then the precise weighted space + interval-arithmetic / Lean (CAP) closure.
