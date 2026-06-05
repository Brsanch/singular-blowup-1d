# Finding: Theta + mu,nu folded in — the convective weighted space IS the symmetry
# quotient; the full frozen shape operator is gapped (complete frozen-sector verdict)

Date 2026-06-05. Code: `validate/theta_munu_gap.py` (PASS), on `src/stability_edge_z.py`.
Completes the frozen (mu=nu=0 base) linear stability verdict for Conjecture 2.4 by
folding in the Theta sector and the mu,nu scaling/normalization modes.

## The clean structural point

Rather than add `mu,nu` with a paper-specific normalization, the convective weighted
space `L^2(e^{-2 gamma z} dz)` (the one giving the spectral gap, see
`findings_spectral_gap.md`) AUTOMATICALLY excludes the scaling symmetry modes:

- **amplitude (`nu`) generator** `dOmega=Omega_bar, dTheta=2 Theta_bar` => in the
  weighted variables `G=(Y-1)^{1/2} dOmega = const`, `h=dTheta = const`. A constant is
  non-normalizable in `L^2(e^{-2 gamma z} dz)` (diverges at the edge `z->-inf`), so it
  is excluded by the weight. It is exactly the marginal `k=0` mode.
- **dilation/edge generator** `X Omega_bar_X ~ (X-1)^{-3/2}` (more singular than the
  profile) -> excluded by the weighted space, and it is the appended `xi` mode at
  `lambda=c_l-1=1` handled separately.

So "asymptotic stability modulo the symmetry group" is realized concretely: the
symmetry modes are precisely the directions the weighted space excludes, and on the
weighted space the full shape operator is gapped. No paper-specific normalization is
needed for the verdict (the explicit-`mu,nu` route should agree).

## Results (validated)

`validate/theta_munu_gap.py`:
1. **Theta sector (`M_hh`)** has the SAME convective gap: weighted-norm decay rate
   `0.500 / 1.000` at `gamma=0.25 / 0.50` (= `2 gamma` exactly); unweighted marginal
   (`-0.000`). Its marginal constant-`Theta` mode is removed by the weight.
2. **Full coupled (`G,h`) operator** `[[M_GG, M_Gh],[0, M_hh]]`: `max Re = -0.000`,
   `#Re>0 = 0`; combined weighted-norm decay `0.497 / 1.000` (= `2 gamma`). Block-
   triangular, so `eig(full) = eig(M_GG) U eig(M_hh)`, both gapped.
3. **Amplitude (`nu`) mode = constant**: weighted norm `~1.1e3` (edge-dominated,
   diverges as `z_min->-inf`) vs a genuine localized shape perturbation `~0.72`
   (finite). The weight excludes the symmetry mode and retains real perturbations.

## Verdict (complete frozen-sector)

The full frozen linearization (Omega-shape + Theta-shape, coupled), in the convective
weighted space, has **no Re>0 point mode and a spectral gap `2 gamma`**; the scaling
symmetry modes (`nu`=constant, dilation/edge=`xi@1`) are exactly the weight-excluded
directions. This is the complete frozen-sector statement of asymptotic stability
modulo symmetries — **numerical evidence consistent with Conjecture 2.4**.

## Honest status / not claimed

- **"Frozen" caveat**: I do not dynamically solve for `mu,nu` via the paper's
  normalization; I argue the weighted space achieves the equivalent symmetry quotient
  (the symmetry modes are non-normalizable there). The explicit-normalization route is
  the referee-facing alternative and should agree; worth doing for a paper.
- **Convective, weighted, numerical.** Not absolute `L^2(dY)` decay (marginal by
  convection); not a proof; ~1.5% Hilbert.
- Remaining: nonlinear shape test on the z-grid; then the precise weighted space and
  the CAP/Lean closure.
