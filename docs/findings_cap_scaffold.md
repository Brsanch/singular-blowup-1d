# Finding: CAP scaffold for chip #3 — the energy/dissipativity route is RULED OUT (real
# non-normality); the CAP must use eigenvalue-enclosure. Obstruction pinned to the
# nonlocal-coupling x weight interaction.

Date 2026-06-06. Code: `validate/cap_scaffold.py` (PASS). First step of the chip-#3 rigor
attempt (certify `sup Re λ < 0` for the linearized stability operator, modulo symmetries).
This scaffold determines WHICH computer-assisted-proof route is viable and pins the leap.

## What was tested

The cleanest CAP for a stability bound is **dissipativity in a weighted inner product**:
if the field of values (numerical range) of the operator `A` in `<u,v>_w` lies in
`Re < 0`, then `sup Re spec(A) < 0` follows, and it reduces to a SYMMETRIC eigenvalue
bound that interval arithmetic can certify. We test this for the Ω-shape operator
`M_GG = transport + nonlocal forcing` in the convective weight `w = e^{-2γz}` (the weight
that gives the spectral gap). Field-of-values abscissa computed via the similarity
transform `Ã = diag(e^{-γz}) A diag(e^{+γz})`, `λ_max((Ã+Ãᵀ)/2)` (numerically clean; no
ill-conditioned generalized eigenproblem).

## Results

- **Transport, analytically (clean):** `<g, −2∂_z g>_w = −2γ‖g‖_w² − [w g²]_{zmin}^{zmax}`.
  Inflow BC `g(zmin)=0` kills the lower boundary term; outflow makes the upper `≤0`. So the
  transport contributes `≤ −2γ‖g‖_w²` (dissipative at rate `2γ`). NB the *discrete*
  symmetrization of `−2∂_z` is boundary-contaminated (FoV abscissa `+33`, not `−2γ`) — the
  one-sided end differences amplified by the edge weight; the transport MUST be handled
  analytically, not by symmetrizing the matrix.
- **Nonlocal forcing `F=½Q_op` (the decisive obstruction):** its weighted FoV abscissa
  `C_force(γ)` is large and **GROWS with edge depth** — `C_force(γ=0.4) = 10.3 → 29.6 → 83.0`
  as `zmin = −10 → −12 → −14`. It does NOT converge to a bound. Reason (real, not a
  discretization artifact): in `<.>_w` the similarity transform multiplies the nonlocal
  forcing entries by `e^{γ(z_j − z_i)}`, which is huge for an edge row (`z_i` small) coupled
  to the bulk — **the exponential weight and the edge↔bulk Hilbert coupling do not commute.**
  So the numerical RANGE in the convective weight is enormous (`Re` up to `+10³⁺`) while the
  SPECTRUM is at `≈ −1e-3`: textbook strong non-normality.

## Determination

**The energy / dissipativity CAP route is RULED OUT in the convective weight.** No `γ`
makes the operator dissipative there (the forcing's weighted numerical range diverges with
edge depth). The transport + boundary signs ARE closed analytically; the obstruction is now
pinned precisely to the **nonlocal-coupling × exponential-weight interaction**.

The CAP must therefore take the **eigenvalue-enclosure route on the spectrum directly** (the
spectrum is benign — `≈ −1e-3`, no `Re>0` mode across two discretizations and the full
nonlinear program; only the field-of-values is bad). Concretely the leap is one of:
1. **Verified eigenpair enclosures** for the low modes (Newton–Kantorovich / Yamamoto a
   posteriori) + a separation bound keeping the rest of the spectrum off `Re=0`, modulo
   `span{ξ@1, μ, ν}`. ⚠️ For a NON-NORMAL operator the a posteriori bound needs resolvent /
   pseudospectral estimates, not the self-adjoint Kato–Temple/Lehmann–Goerisch directly.
2. **Verified resolvent bound** `‖(L−λ)^{-1}‖ ≤ M` on `Re λ ≥ 0` (modulo symmetries) — the
   non-normality-robust route. ⚠️ must be in the WEIGHTED space (the essential spectrum sits
   on `Re=0` unweighted; the gap is only after the weight), so it inherits the same
   `e^{γz}` conditioning that must be controlled with interval arithmetic, not floats.
3. **A non-exponential metric** that bounds the nonlocal term while keeping transport
   dissipative — open; the convective weight provably does not (this scaffold).
Plus, in every case: a rigorous Hilbert-kernel enclosure (the ~1.5% quadrature → verified)
and a truncation/tail bound for the mapped infinite domain; and folding in `M_Gh/M_hh`
(Θ couples via a bounded `e^{-z/2}∂_z` term — same treatment).

## Honest status

This is a CAP *scaffold*, not a CAP. It rigorously (analytically) closes the transport +
boundary part and numerically determines that the energy route fails and why, pinning the
remaining leap to the nonlocal-coupling × weight interaction. The verified-enclosure
execution (routes 1–3) is the human-leap and needs interval-arithmetic tooling (the machinery
De Huang's group owns; scoop-risky). Numerical scaffold; not a proof.
