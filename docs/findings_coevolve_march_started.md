# Finding: co-evolving march STARTED — dynamic-normalization machinery validated
# (profile is an exact fixed point); the fixed-Y wall quantified → moving frame required

Date 2026-06-06. Code: `validate/coevolve_march.py` (PASS). First step of THE NEXT
CHIP #2 (co-evolving-base dynamic-rescaling march with explicit normalization). This
is a START, not the completed relaxation result: it builds + validates the
dynamic-normalization half, and pinpoints/quantifies the obstruction that dictates
the frame for the rest.

## What this chip is (vs the frozen-base march already done)

The frozen-base nonlinear march (`findings_nonlinear_shape_march.md`) evolved the
perturbation at FIXED `c_l=2, c_omega=-1` — that is exactly the full δΩ-dynamics with
the exponents held. The co-evolving chip lets `c_l, c_omega` adjust DYNAMICALLY by the
normalization and asks: does a perturbed profile relax to `(Ω̄, c_l=2, c_omega=-1)`
with no exponent drift? Per `findings_dynamic_rescaling.md` the profile is a SADDLE
whose single unstable direction is the amplitude; the relaxation requires a
normalization that pins the gauge (projects out that direction).

## Built + validated: the dynamic normalization (split form)

Two conditions, derived so that both collapse to the exact profile values:
- `c_l = -U(1)` (edge-stagnation: `ā(1)=U(1)+c_l=0`, pins the edge at Y=1) → `c_l=2`.
- `c_omega = [ā·Ω_Y - Θ_Y]/Ω` at a bulk reference `Y_ref` (amplitude pin) → `c_omega=-1`.

Computed in the VERIFIED split (closed-form `Ūbar=-2` on Y>1 + bounded perturbation
Hilbert on `G=(Y-1)^{1/2}δΩ`), NOT by Hilberting the raw singular profile (which is
~e^{|zmin|/2} at the edge and under-resolves). The only new numeric piece is
`dU(1)=∫_0^1 H(δΩ)dY'` (the perturbation velocity reaching the edge from the support,
a regular `Y'<1<Y` integral):

- **GATE 0a**: `dU(1)` functional `b1@G` vs direct fine quadrature — rel.err **7.6e-8**.
- **GATE 0b**: on the exact profile (`G=0`): `dU(1)=0`, `c_l=2.00000000`,
  `c_omega=-1.00000000` (machine precision), `ā(edge)≈1e-5`. **The profile is an exact
  fixed point of the normalization.** (An earlier raw-Hilbert attempt gave 2.24 / -1.11
  and blew up even from the profile — confirmed to be an edge under-resolution artifact,
  now eliminated.)

## The wall, quantified: fixed-Y frame is invalid → moving frame mandatory

When `c_l ≠ 2` (any perturbation gives `δc_l = -dU(1) ≠ 0`), eq.(2.4) acquires the
gauge-correction forcing `-δc_l·X·Ω_X`. Its base part `-δc_l·X·Ω̄_X`, in the bounded
G-variable, is `(δc_l/2)(1+e^z)e^{-z} ~ (δc_l/2)e^{-z}` near the edge. Measured edge
magnitude (bulk bump, `δc_l=+0.044` fixed):

| zmin | N | δc_l | max\|gauge term\| | at edge node |
|---|---|---|---|---|
| -10 | 1000 | +0.0438 | 4.82e2 | +4.82e2 |
| -12 | 1400 | +0.0438 | 3.56e3 | +3.56e3 |
| -14 | 1800 | +0.0438 | 2.63e4 | +2.63e4 |
| -16 | 2200 | +0.0438 | 1.95e5 | +1.95e5 |

The term **diverges as `e^{-z}`** (each `Δzmin=-2` multiplies by `e²≈7.4`: 482→3563→
26330→194500). Root cause: `X·Ω̄_X ~ (Y-1)^{-3/2}` is MORE singular than the profile
`(Y-1)^{-1/2}` — it is precisely the dilation/edge symmetry generator that the linear
verdict's weighted space EXCLUDES (`findings_theta_munu_gap.md`). A co-evolving march
therefore CANNOT run in the fixed-Y frame: the `c_l` (dilation) gauge produces a
spurious infinite edge velocity (the same mechanism behind the "δ-divergent edge
artifact" recorded in the dead-ends).

## Verdict / status

- **Validated**: the dynamic-normalization machinery (edge-stagnation `c_l` + amplitude-
  pin `c_omega`) is correct and makes the profile an exact fixed point.
- **Pinpointed + quantified**: the fixed-Y frame is invalid for the co-evolving march
  because the `c_l` gauge couples to the edge-singular `(Y-1)^{-3/2}` mode. **The moving
  edge-fixed frame (ξ-ODE) is mandatory** — the same formulation the verified LINEAR
  operator `M_GG` already uses (`findings_edge_condition_derivation.md`).
- This is a START. NOT a relaxation result; not a proof.

## Precise next step (specified)

Build the co-evolving march in the MOVING frame:
`co-evolving march = nonlinear_shape_march (L = M_GG,M_Gh,M_hh + nonlinear N) +
ξ-ODE  ξ̇ = (c_l-1)ξ + V_g(1) + μ  (V_g(1)=∫_0^1 H(g)dY = the `b1@G` built here) +
explicit μ,ν normalization (amplitude/scaling).` In this frame the dilation is absorbed
by ξ (no singular term), so `c_l, c_omega` can co-evolve. Confirm the perturbed profile
relaxes to `(Ω̄, c_l=2, c_omega=-1)` with no exponent drift under refinement; HARD STOP
if exponents drift. ⚠️ the paper's explicit μ,ν normalization conditions should be
matched (referee-facing) and confirmed to agree with the weighted-space quotient.
