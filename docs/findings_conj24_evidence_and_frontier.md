# Conjecture 2.4 — consolidated numerical evidence (frozen sector) + the precise rigor frontier

Date 2026-06-06. Capstone for the stability program on the singular Hou–Luo profile
(arXiv:2604.01868 Thm 2.3, verified real). The numerical case is now comprehensive; chip #3
(rigor / CAP) is the human-leap. This doc assembles the evidence and states the exact
obstruction at maximum resolution.

## The target

Conjecture 2.4 (asymptotic stability of the singular self-similar profile
`Ω̄=1{X>1}/√(X-1)`, `Θ̄=(π/2)1{X>1}`, `c_l=2`, `c_ω=-1`), in the convective weighted sense
standard for transport-dominated self-similar operators (Chen–Hou / De Huang lineage) —
NOT absolute `L²(dY)` decay. The paper is numerical-only and publishes no stability spectra;
this is the hard-to-scoop contribution.

## Evidence chain (each item independently validated, with convergence + resolution checks)

| # | result | gate |
|---|---|---|
| 1 | profile certified (closed form, `U(1)=-2`, residual 2e-16) | `verify_singular_steady.py` |
| 2 | linear spectrum (frozen): NO `Re>0` point mode; rightmost `≈-0.001` = marginal transport continuum; two independent discretizations; edge-resolved, no box modes | `stability_edge_z_spectrum.py` |
| 3 | convective spectral GAP `2γ` in `L²(e^{-2γz}dz)`; symmetry modes (ν=const, dilation/edge=ξ@1) are exactly the weight-excluded directions | `spectral_gap.py`, `theta_munu_gap.py` |
| 4 | operator = correct linearization of the trusted eq.(2.4) flow (rel.err ∝ ε→0) | `nonlinear_linearization_check.py` |
| 5 | no transient growth for physical perturbations (`G=1`); non-normal worst case is a grid-scale edge artifact | `transient_growth.py` |
| 6 | nonlinear FINITE-amplitude (frozen base): convective decay to `A=4`, no transient growth | `nonlinear_shape_march.py` |
| 7 | λ=1 pure-ξ symmetry mode confirmed in the coupled (g,h,ξ) operator | `coevolve_moving_frame.py` |
| 8 | fully self-consistent CO-EVOLVING march relaxes: both `c_l,c_ω` fed back → `(2,-1)`, no drift, 6-digit resolution-robust | `coevolve_projected.py` |
| 9 | coupled (Ω,Θ)-shape relaxes; stress test (amplitude≤2; bulk/near-edge/far/high-k/two-bump/Θ-dominant) finds NO counterexample | `coevolve_projected.py`, `coevolve_stress.py` |
| 10 | normalization matches the paper's eq (3.2) (`c_l=-U(1)` identical; `U_x(0)→-1`; exponents `(2,-1)`); both anchor papers verified real | `coevolve_paper_match.py` |

**Net:** linear stability (no unstable mode + convective gap + no transient growth) AND
nonlinear stability (finite-amplitude + fully co-evolving + coupled + stress-tested) AND
formulation correctness (correct linearization + paper-matched normalization). Coherent,
multi-method, no surprises.

## The precise rigor frontier (chip #3 — the obstruction at maximum resolution)

**Exact statement to certify:** the global linearized operator `L` (shape eqs §4 of
`findings_edge_condition_derivation.md`, with the NONLOCAL `∫_1^Y H(g)` coupling, on a
mapped infinite `Y>1` domain), restricted to the convective weighted space
`L²(e^{-2γz}dz)` and quotiented by the symmetry subspace `span{ξ@λ=1, μ, ν}`, has
**`sup Re λ < 0`** (in fact the essential spectrum is at `Re=-2γ`; the claim is that NO
point mode crosses into `Re≥0`).

**Why numerics cannot close it:** finite discretization + a ~1.5% line-Hilbert gives the
spectral abscissa to ~1e-3, not a verified sign. The marginal `≈-0.001` rightmost eigenvalue
is the discretized transport continuum, not a certified bound. Closing it needs a
**computer-assisted proof** (interval arithmetic on the operator + rigorous Hilbert-kernel
and edge enclosures), the machinery De Huang's group owns — hence scoop-risky.

**What's encouraging (entry points):**
- The local edge indicial analysis FORCES any `Re λ>0` mode to be non-singular at the edge
  (drop the homogeneous `(Y-1)^{-(λ+1)/2}` branch; Robin BC) — it cannot hide in the edge.
- No transient growth for physical perturbations (`G=1`) bounds the non-normal amplification
  that could otherwise seed a nonlinear instability.
- The only residual risk is the **nonlocal `∫H(g)` coupling** carrying a point mode into
  `Re≥0` — but items 2/3/8/9 see none, across two discretizations, finite amplitude, and a
  stress sweep.

**The leap (handoff):** turn "no `Re λ>0` numerically" into a CAP-certified
`sup Re λ < 0` (modulo symmetries) on the nonlocal operator in the weighted space — or find
the analytic structure that rules out a nonlocal-coupling unstable mode. This is the
non-interpolative step; the numerics have made it well-posed and sharply localized (it is
purely about the nonlocal coupling, the edge and transient-growth risks being closed).

## Honest status

Numerical evidence consistent with Conjecture 2.4 across linear + nonlinear + co-evolving +
paper-matched checks; convective/weighted; frozen base resolved (co-evolving done); NOT a
proof. The rigor closure (chip #3) is open and scoop-risky.
