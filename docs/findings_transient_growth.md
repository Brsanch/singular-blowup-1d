# Finding: transient growth / non-normality is MODEST (supports robustness) but the
# precise bound is edge-conditioning-limited (an honest open refinement)

Date 2026-06-05. Code: `validate/transient_growth.py`. The question beyond the
eigenvalues: the operator is strongly NON-NORMAL (transport-dominated), so stable
eigenvalues + a gap do not by themselves preclude large TRANSIENT amplification
`||e^{M tau}||`, which is what can trigger NONLINEAR instability. Is it modest or large?

## Conditioning had to be fixed first (a real trap)

The naive computation gave absurd, resolution-dependent amplification (millions, no
turnover, numerical abscissa ~+8500) -- a discretization artifact: the forcing
`(1/2) J/(Y-1)` is analytically bounded (`-> (1/2) H(g)(1)`) but discretely it is
`small/small` with `1/(Y-1)=e^{-z}` huge at the edge. Recasting it as the manifestly
bounded exponential-averaging convolution
`(1/2) J/(Y-1) = (1/2) int_{zmin}^z H(g) e^{z'-z} dz'` (weights `<= 1`) makes `||M||`
RESOLUTION-INDEPENDENT (~220-260), so amplification becomes trustworthy.

## Result (physical L^2(dY) norm)

| N | zmin | G_max = sup_tau ||e^{M tau}|| | ||e^{12 M}|| |
|---|---|---|---|
| 600 | -8 | 1.76 (tau=0.5) | 1.00 |
| 800 | -10 | 2.33 (tau=0.25) | 1.00 |
| 1000 | -12 | 3.08 (tau=0.25) | 1.00 |

- The transient amplification is **MODEST (G_max ~ 2-3)** and **decays back toward
  ~1** (bounded transient, NO blow-up) -- supportive of robust linear stability and
  hence nonlinear stability.
- **Caveat (honest):** `G_max` creeps up with edge depth (`1.8 -> 2.3 -> 3.1` as
  `zmin -> -12`), which cannot be cleanly separated from a residual `e^{-z/2}` edge
  factor in the Hilbert map `H(g)`. So the PRECISE non-normality bound is an OPEN
  refinement (needs a fully edge-conditioned operator). The weighted worst-case
  operator norm is also conditioning-contaminated (far-field `e^{gamma z}` weight),
  so it is not reported; the convective DECAY of specific perturbations is the
  reliable weighted-norm statement (`findings_spectral_gap.md`).

## Net (calibrated, not overclaimed)

- **Solid:** eigenvalue stability (no Re>0 point mode, two grids), the convective
  spectral gap (specific-perturbation semigroup decay), the operator = correct
  linearization of the trusted nonlinear flow.
- **Suggestive but not settled:** the transient-growth / non-normality bound is
  modest (~2-3, decaying) in the best computations, which SUPPORTS nonlinear
  stability, but the resolution-sensitivity means it is not a clean robustness proof.
  A fully edge-conditioned operator (bounded `H(g)` map; e.g. an analytic edge
  treatment or a `(Y-1)^{1/2}`-symmetrized variable) would settle it.

This is the current frontier of the numerical stability study: the linear verdict is
robust; the non-normal transient-growth robustness is the remaining numerical
sub-question (and the natural precursor to a rigorous nonlinear-stability statement).
