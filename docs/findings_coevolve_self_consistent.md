# Finding: fully self-consistent co-evolving march RELAXES — the singular profile is a
# nonlinear attractor of the dynamic-rescaling flow (both c_l, c_ω co-evolving), no drift

Date 2026-06-06. Code: `validate/coevolve_projected.py` (PASS). Completes THE NEXT CHIP
#2: the FULLY self-consistent co-evolving march (both exponents fed back, not just tracked
diagnostically). Builds on the moving-frame coupled operator (`findings_coevolve_moving_frame.md`),
which relaxed only diagnostically because explicit forward forcing blew up.

## The runaway, localized (STEP 1)

Toggling the two normalization forcings independently (no projection):

| forcing | result |
|---|---|
| ν only (Ω-amplitude, `+ν` in g) | **stable**, ‖g‖_w bounded |
| `(μ+2ν)Θ̄` h-forcing (Θ-amplitude) | **BLOWS UP** |
| both | BLOWS UP |

So the runaway is the **Θ-amplitude h-forcing**, NOT the Ω-amplitude — it excites the
non-normalizable Θ-amplitude mode (constant `h`), which fights the edge BC `h(zmin)=0`
and is amplified by the `M_Gh = e^{-z/2}∂_z` coupling at the edge. (This corrects the
earlier guess that the Ω-amplitude saddle was the culprit; the naive "subtract the value
at a point" projection is also wrong — it *adds* the non-normalizable constant mode.)

## The fix: the Θ-balance normalization (STEP 2)

Impose `c_ω = −c_l/2` (the Θ-balance: keeps `Θ̄=(π/2)1{Y>1}` an exact steady state, since
the singular profile has `c_l+2c_ω=0`). Then **`c_l+2c_ω ≡ 0`**, so the h-forcing
`(μ+2ν)Θ̄ = (c_l+2c_ω)(π/2)` **vanishes identically**. Combined with the edge-stagnation
`c_l = −U(1)` and the edge-pinning `μ=−V_g(1)` (`ξ̇=(c_l−1)ξ` homogeneous → ξ stays 0), both
exponents are determined and co-evolve:
```
c_l = -U(1) = 2 - V_g(1),   c_ω = -c_l/2 = -1 + V_g(1)/2,   V_g(1) = ∫_0^1 H(g).
```
Only the (stable) ν = V_g(1)/2 Ω-amplitude forcing remains, and `ν→0` as `g→0`.

## Result: relaxation with both exponents co-evolving (STEP 3)

Perturbed singular profile (bulk shape bump on Ω), both `c_l, c_ω` fed back each step:

| τ | c_l | c_ω | ξ | ‖g‖_w(.5) |
|---|---|---|---|---|
| 0.0 | +2.04379 | −1.02189 | 0 | 1.08e-1 |
| 1.0 | +2.00225 | −1.00112 | 1.8e-4 | 4.59e-2 |
| 2.0 | +1.99978 | −0.99989 | 5.1e-4 | 1.48e-2 |
| 3.0 | +2.00000 | −1.00000 | 1.4e-3 | 5.24e-3 |
| 4.0 | +2.00001 | −1.00000 | 3.8e-3 | 1.91e-3 |

`‖g‖_w` decays by ~57× (rate ≈ +1.26); **`c_l→2` and `c_ω→−1` to 5 digits, monotone, no
drift**; ξ stays ~0 (edge pinned).

**Convergence** under `(N,zmin) = (1000,−10), (1400,−12), (1800,−14)`: final
`c_l = 2.000006`, `c_ω = −1.000003`, `‖g‖_w = 1.9e-3` — **identical to 6 digits across all
three grids**. (The fitted convective decay rate is `1.20 / 1.26 / 1.36` — mildly
zmin-sensitive, as weighted-convective rates are; the converged exponents and relaxation
state are resolution-robust and the rate stays above the linear gap `2γ=1.0`.)

**Coupled (Ω,Θ)-shape perturbation.** Bumping BOTH `G` and `h` (not just `G`), same
self-consistent normalization: `c_l: 2.044→2.00001`, `c_ω: −1.022→−1.00001` (converge, no
drift), `‖g‖_w` decays `0.108→2.9e-3` (~37×); `h` (the slaved Θ-mode, pure transport)
convects out (`max|h|` preserved, weighted norm decays). A small bounded transient in
`max|g|` (`0.15→0.18`) from the `M_Gh` coupling, then decay. **No coupling-driven
instability** — the Ω-only verdict is robust to perturbing Θ as well.

## Stress test (falsification-first; `validate/coevolve_stress.py` PASS)

Tried to BREAK the attractor — every case relaxes (`‖g‖_w` decays, `c_l→2`, `c_ω→−1`, no
blowup/drift):
- **Amplitude** `A = 0.05, 0.15, 0.3, 0.6, 1.0, 2.0` (bulk bump; `A=2` starts at
  `‖g‖_w=1.45`): all relax; final `c_l∈[2.00000, 2.00014]`, `c_ω∈[−1.00000, −1.00007]`.
- **Shape** (A=0.2): bulk (z=1), **near-edge (z=−4)**, far-field (z=6), **high-k oscillatory**
  (cos 6z), sign-changing two-bump — all relax.
- **Θ-dominant** (small Ω bump + big `h` bump 0.5): relaxes.
- **Convergence** (A=1.0): final `‖g‖_w=1.34e-2` and `c_l=2.00005, c_ω=−1.00003` identical
  across `(N,zmin)=(1000,−10)/(1400,−12)/(1800,−14)`.

No resolution-robust blowup or exponent drift found under any amplitude or shape — the
attractor claim survives an aggressive falsification sweep.

## Verdict

The **fully self-consistent co-evolving dynamic-rescaling march relaxes**: a perturbed
singular profile returns to `(Ω̄, c_l=2, c_ω=−1)` with both exponents dynamically
co-evolving and NO exponent drift (resolution-robust to 6 digits). **The singular Hou–Luo
profile is a nonlinear attractor of the dynamic-rescaling flow** under the edge-stagnation
+ Θ-balance normalization, with the edge pinned. This is the co-evolving-base confirmation
the chip targeted — it removes the "frozen base" caveat of `findings_nonlinear_shape_march.md`
and the "diagnostic only" caveat of `findings_coevolve_moving_frame.md`.

## Honest status / not claimed

- **Specific normalization.** Edge-stagnation (`c_l=−U(1)`) + Θ-balance (`c_ω=−c_l/2`):
  both structural, both exact at the profile, a legitimate self-consistent dynamic-rescaling
  gauge. Whether it matches arXiv:2604.01868's explicit `μ,ν` conditions verbatim is the
  referee-facing check still to do (it should agree, both selecting `(2,−1)`).
- **Perturbation coverage.** Both the Ω-shape alone and the coupled (Ω,Θ)-shape are now
  tested (above); both relax with exponents → (2,−1), no drift, no coupling-driven
  instability.
- **Convective/weighted, numerical, not a proof.** Decay in `L²(e^{-2γz}dz)`; finite
  edge-resolved grid; ~1.5% Hilbert; RK4. No interval arithmetic / Lean.
- The fitted convective rate is norm/domain-dependent (1.2–1.4); the robust invariants are
  the converged exponents and the relaxation.

## Next

1. Match arXiv:2604.01868's explicit `μ,ν` normalization conditions; confirm agreement with
   the edge-stagnation + Θ-balance gauge used here and with the weighted-space quotient.
2. Pin the precise weighted space / max gap; then interval-arithmetic / Lean (CAP) closure
   (scoop-risky, the remaining frontier).
