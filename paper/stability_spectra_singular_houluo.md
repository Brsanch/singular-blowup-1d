# Stability spectra of the singular self-similar profile of the 1D Hou–Luo model

> **SUPERSEDED (2026-07-22).** Publish-ready package is
> `paper/stability_spectra_singular_houluo.tex` → **`.pdf`** (ready-to-arXiv).
> Figures: `paper/scripts/fig_*.py`. Scoop log:
> `docs/scoop_recheck_2026_07_22_stability_spectra.md`.
> Keep this file as historical draft only.

**Draft, 2026-06-06 (historical).** Companion to the numerical-existence papers of Chen–De Huang–Li
(arXiv:2604.01868) and Huang–Tong–Wang (arXiv:2603.25104), which construct the singular
(unbounded) self-similar blowup profiles but provide **no stability analysis**. This draft
reports the first **stability spectra** for the singular Hou–Luo self-similar profile —
linear and nonlinear, numerical — and the analytical estimate to which a rigorous proof
reduces. All claims are reproducible by the `validate/` scripts named inline.

---

## Abstract

For the 1D Hou–Luo model `ω_t + u ω_x = θ_x`, `θ_t + u θ_x = 0`, `u_x = H[ω]`, the singular
self-similar profile `Ω̄(X)=1{X>1}/√(X−1)`, `Θ̄=(π/2)1{X>1}` (`c_l=2`, `c_ω=−1`; Thm 2.3 of
arXiv:2604.01868) is shown numerically to be **asymptotically stable modulo its symmetry
group**, in the convective weighted sense standard for transport-dominated self-similar
operators. The linearized operator, desingularized in the edge-resolving coordinate
`z=ln(Y−1)` and in the moving (edge-fixed) frame, has **no eigenvalue with `Re λ>0`** (two
independent discretizations), a **convective spectral gap `2γ`** in `L²(e^{−2γz}dz)`, and the
symmetry modes (time-translation `ξ` at `λ=1`; scaling/amplitude `µ,ν`) are exactly the
weight-excluded directions. **No transient growth** occurs for physical perturbations, and in
the edge-conditioned energy norm the operator is a strict contraction (numerical abscissa
`<0`, resolution-independent). Nonlinearly, a fully self-consistent dynamic-rescaling march
(both rescaling exponents co-evolving by the paper's normalization) **relaxes** a perturbed
profile back to `(Ω̄, c_l=2, c_ω=−1)` with no exponent drift, robust across amplitudes (to
`4×` the profile scale) and perturbation shapes. A weighted-energy finite-codimension
estimate is set up; its key coercivity inequality is verified numerically to hold in the
bulk, reducing a rigorous proof to a matched edge-layer asymptotic plus a higher-order
closure.

## 1. Setup

**Model and dynamic rescaling.** In self-similar variables (eq. (2.4) of arXiv:2604.01868),
```
Ω_τ + (U + c_l X)Ω_X = c_ω Ω + Θ_X,   Θ_τ + (U + c_l X)Θ_X = (c_l+2c_ω)Θ,   U_X = H(Ω), U(0)=0.
```
**Certified singular steady state** (Thm 2.3, re-verified here, `validate/verify_singular_steady.py`):
`Ω̄=1{X>1}/√(X−1)`, `Θ̄=(π/2)1{X>1}`, `Ū(X)=2√(1−X)−2` (`X<1`), `−2` (`X>1`); `c_l=2`, `c_ω=−1`;
steady residual `2e-16`, `U(1)=−2` exactly. The gCLM/De Gregorio engine is validated against
the CLM closed form to `2e-10` (`validate/validate_clm.py`); the physical-space Hou–Luo
solver reproduces `c_ω=−1` to `<0.4%` (`findings_houluo_physical.md`).

**Linearized operator (edge-resolved, moving frame).** Perturbation `Ω=Ω̄+g`, `Θ=Θ̄+h`,
edge displacement `ξ`. In `z=ln(Y−1)` with the bounded variable `G=(Y−1)^{1/2}g`:
```
M_GG = −2∂_z + ½Q_op,   Q_op G = δU/(Y−1),  δU=∫_1^Y H(g);   M_hh = −2∂_z;   M_Gh = e^{−z/2}∂_z;
ξ̇ = (c_l−1)ξ + V_g(1) + µ,   V_g(1)=∫_0^1 H(g).
```
`M_GG = −2∂_z + ½Q_op` verified to `1e-11`; the `c_ω=−1` source is absorbed into the clean
`−2∂_z` transport by the `(Y−1)^{1/2}` conjugation. The edge condition is Robin (drop the
homogeneous `(Y−1)^{−(λ+1)/2}` branch; `findings_edge_condition_derivation.md`), and the
far-field is outflow (non-periodic); the moving frame removes the spurious infinite edge
velocity of the fixed-frame linearization.

## 2. Linear stability spectrum

`validate/stability_edge_z_spectrum.py`, `spectral_gap.py`, `theta_munu_gap.py`,
`transient_growth{,_edge_conditioned}.py`.

- **No unstable point mode.** `max Re eig(M_GG) ≈ −1e-3` (the marginal transport continuum);
  `#{Re>0}=0`, in two independent discretizations (equispaced and edge-resolved `z`), edge-
  resolved with no spurious near-edge cloud and no box modes.
- **Convective spectral gap.** The transport `−2∂_z` is convective: marginal in `L²(dY)`,
  but in `L²(e^{−2γz}dz)` the essential spectrum is at `Re=−2γ`. Semigroup decay rates
  (`spectral_gap.py`): `+0.054 / +0.552 / +1.057 / +1.675` at `γ=0 / 0.25 / 0.5 / 0.75`,
  i.e. `≈ 2γ`, linear in `γ` (no internal point mode caps the gap).
- **Symmetry modes are weight-excluded.** The amplitude generator (`ν`) is `G=const`
  (non-normalizable in the weight); the dilation/edge generator (`µ`) is `~(X−1)^{−3/2}`
  (excluded); the time-translation mode is the pure-`ξ` eigenvalue `λ=c_l−1=1` (confirmed
  pure-`ξ`). The full coupled `(G,h)` operator is gapped with these projected out.
- **No transient growth (settled).** Physical (smooth) perturbations have amplification
  `G=1` at every resolution. In the edge-conditioned energy norm (bounded-at-edge anisotropic
  weight, transport via exact IBP) the numerical abscissa is `µ_an<0` (`−0.10` at `z_cut=−8`,
  `−0.41` at `−6`) and **resolution-independent** ⇒ `‖e^{Mτ}‖_W ≤ 1` (strict contraction).
  The earlier worst-case creep (`3→4` with edge depth) is the discretized `(Y−1)^{−1/2}` edge
  spike, excluded by the Robin condition.

## 3. Nonlinear stability

`validate/nonlinear_shape_march.py`, `coevolve_*.py`.

- **Finite-amplitude (frozen base).** A finite shape perturbation decays in the convective
  weighted norm for amplitudes to `A=4` (`4×` profile scale), with `peak ‖·‖_w/‖·‖_0 = 1.000`
  (monotone, no transient growth) and resolution-robust rate.
- **Fully self-consistent co-evolving relaxation.** Letting both rescaling exponents adjust by
  the normalization (edge-stagnation `c_l=−U(1)` + Θ-balance `c_ω=−c_l/2`, equivalent at and
  near the profile to the paper's eq. (3.2)), a perturbed profile relaxes to
  `(Ω̄, c_l=2, c_ω=−1)`: exponents converge to 6 digits with **no drift**, identical across
  `(N,zmin)` grids; edge pinned. Coupled `(Ω,Θ)` perturbations also relax (no coupling-driven
  instability), and a stress sweep (amplitude `≤2`; bulk/near-edge/far/high-k/two-bump/
  Θ-dominant shapes) finds **no counterexample**.

## 4. Toward a rigorous proof

`docs/rigor_energy_estimate_setup.md`, `validate/energy_estimate_search.py`. A weighted-energy
finite-codimension estimate (à la Chen, arXiv:2605.15130 — anisotropic weights + integration
by parts exploiting the equation twice) reduces stability to the coercivity
`⟨(G,h), 𝓛(G,h)⟩_{H¹_W} ≤ −δ‖(G,h)‖²` for `G,h ⊥ span{ξ,µ,ν}`. The transport contributes
`−∫φ'(z)W G² dz` (exact IBP; inflow boundary vanishes by Robin, outflow `≤0`); the key
inequality is that the nonlocal forcing and the `M_Gh` coupling are dominated by this
dissipation. With a bounded-at-edge anisotropic weight, the **bulk energy abscissa is
negative** (`−0.02 → −0.41` excluding a thin edge layer where `G→0`) — i.e. the coercivity
holds away from the edge. A rigorous proof reduces to: (i) a matched edge-layer asymptotic;
(ii) the `H¹` closure of the `M_Gh` coupling; (iii) the nonlinear step.

## 5. Honest status / limitations

- **Numerical, not a proof.** All stability claims are numerical (finite edge-resolved
  discretization, `~1.5%` line-Hilbert quadrature). §4's coercivity is verified numerically,
  not proven analytically; the rigorous estimate is open (the analytic leap).
- **Convective / weighted.** Stability is asymptotic in `L²(e^{−2γz}dz)` modulo symmetries —
  the standard sense for transport-dominated self-similar operators, **not** absolute `L²(dY)`
  decay (which is marginal by convection).
- **Scope.** 1D Hou–Luo singular profile; the 2D Boussinesq analogue and the gCLM `a`-family
  are out of scope here.

## References

- B. Chen, D. Huang, X. Li, *Novel self-similar finite-time blowups … 1D Hou–Luo and 2D
  Boussinesq: a numerical investigation*, arXiv:2604.01868 (2026).
- D. Huang, J. Tong, X. Wang, *Self-similar finite-time blowups … generalized
  Constantin–Lax–Majda model*, arXiv:2603.25104 (2026).
- J. Chen, *Asymptotically self-similar blowup for 3D incompressible Euler with `C^{1,1/3−}`
  velocity I, II*, arXiv:2605.15149, 2605.15130 (2026) — finite-codimension stability template.

*Reproducibility: every numerical claim is produced by a `validate/*.py` gate (pure Python,
1D, cheap); see `HANDOFF.md` for the full state-of-play table.*
