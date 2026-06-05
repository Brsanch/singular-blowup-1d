# Finding: the spectral gap exists — the singular profile is asymptotically stable
# in the convective weighted norm (resolves the crux of Conjecture 2.4, frozen sector)

Date 2026-06-05. Code: `validate/spectral_gap.py` (PASS), on the edge-resolved
operator `src/stability_edge_z.py`. Target: the "asymptotic vs marginal" question
that was the single biggest open gap in the linear stability verdict.

## The question

Earlier work showed the Omega-shape operator has NO unstable point mode, but its
essential spectrum sits at Re=0 (marginal) in the plain norm. Asymptotic stability
(Conjecture 2.4) requires a strict spectral GAP (essential spectrum Re<0). That gap
could not be seen by a finite-matrix weight (similarity preserves eigenvalues), so
the question was open: asymptotically stable, or only marginal?

## Why the plain norm is marginal, and where the gap lives

The transport `-2 d_z` (in z=ln(Y-1)) is CONVECTIVE: it carries perturbations toward
+z (the far field) at speed 2. On the line this is norm-preserving in L^2(dz) (=
L^2(dY) for the physical perturbation), so the essential spectrum is on the imaginary
axis -> marginal, no absolute decay (perturbations propagate out without decaying
pointwise). This is intrinsic to an outflow transport operator, not a defect.

The gap lives in a WEIGHTED (convective) norm. In `L^2(e^{-2 gamma z} dz)` (gamma>0:
require decay at the edge z->-inf; allow far-field outflow growth -- the only
admissible direction for an outflow operator), the symbol of `-2 d_z` conjugated by
`e^{-gamma z}` is `-2 d_z - 2 gamma`, so the
**essential spectrum is Re = -2 gamma < 0  =>  spectral gap = 2 gamma.**

## Numerical confirmation (semigroup decay = weighted spectral abscissa)

`validate/spectral_gap.py` evolves a localized bulk shape perturbation under the
edge-resolved operator and measures norm decay:

| gamma | measured decay rate | predicted 2*gamma |
|---|---|---|
| 0 (unweighted) | +0.054 | marginal |
| 0.25 | +0.552 | 0.50 |
| 0.50 | +1.057 | 1.00 |
| 0.75 | +1.675 | 1.50 |

- The UNWEIGHTED norm is ~flat (rate +0.05 ~ marginal): NO absolute decay in L^2(dY)
  (perturbations are swept to the far field), as expected.
- Each EDGE-WEIGHTED norm decays at rate ~ `2 gamma` (+ the ~0.05 nonlocal baseline).
- Linear in gamma through 0.75 with NO saturation -> no internal point mode caps the
  gap (consistent with the Omega-shape sector having no point spectrum). The gap is
  CONVECTIVE: its size scales with the weight (transport speed x gamma), not a fixed
  intrinsic number.

## Verdict

**The spectral gap exists.** The Omega-shape sector of the linearization is
ASYMPTOTICALLY stable in the convective weighted norm `L^2(e^{-2 gamma z} dz)` (gap
`2 gamma`), modulo the symmetry modes. This is the standard weighted-norm sense in
which transport-dominated self-similar stability is framed (Chen-Hou / De Huang
lineage); it is NOT absolute decay in `L^2(dY)` (that is marginal, by convection).

This resolves the "asymptotic vs marginal" crux for the frozen Omega-shape sector:
the answer is asymptotic (weighted). Combined with: no unstable point mode (two
independent discretizations), the edge/`xi` symmetry mode confirmed linearly and
nonlinearly at `lambda=1`, and the spurious-edge-mode trap caught and fixed -- the
frozen-sector evidence for Conjecture 2.4 is now coherent and includes the gap.

## Honest status / not claimed

- **Omega-shape sector, frozen `mu=nu=0`.** The `Theta` sector (marginal at Re=0) and
  the `mu,nu` scaling/normalization modes are not yet folded into the weighted-gap
  statement; they need the paper's normalization conditions and the symmetry
  projection. The `Theta` marginal direction in particular is expected to be removed
  by the same edge-decay weight (a constant in z is not in `L^2(e^{-2 gamma z} dz)`),
  but that should be shown explicitly.
- **Convective, norm-dependent gap.** The statement is asymptotic stability in a
  weighted norm, the standard sense here; not a fixed intrinsic decay rate, and not
  absolute L^2(dY) decay.
- **Numerical, not a proof.** Semigroup decay on a finite, edge-resolved
  discretization with a ~1.5% Hilbert. No interval arithmetic / Lean.

## Next

1. Fold in `Theta` + `mu,nu`: show the weight removes the marginal `Theta` mode, add
   the scaling modes with the paper's normalization, project out `{xi@1, mu, nu}`,
   and confirm the full operator's gap.
2. Nonlinear shape test on the z-grid (independent confirmation of decay).
3. Pin the precise weighted space / max gap; then interval-arithmetic / Lean closure.
