# singular-blowup-1d

Self-similar **finite-time blowup** of 1D nonlocal active-scalar models
(Constantin–Lax–Majda / De Gregorio / generalized-CLM / Hou–Luo), with a focus
on the **singular-profile** scenarios discovered in 2026 and the rigorous /
computer-assisted certification of their profiles and **linear stability**.

## The bet (why this, and why now)

Years of work on Navier–Stokes / SQG / Euler **global regularity** hit a hard
math wall: every local geometric-depletion route collapses to Type-I /
KNSS–Liouville rigidity — a statement equivalent to the conclusion, not a
measurable gap. This repo **inverts that wall**. Instead of trying to *rule out*
blowup on the true equation, we *certify* blowup on toy models where blowup is
the answer. These are 1D, so extreme spatial resolution is cheap, and the
profile ODEs are the cleanest computer-assisted-proof (CAP) / Lean targets in
the whole program. None of the listed dead ends apply: a 1D model is not 3D
Euler regularity, not Type-I rigidity, not undecidable.

## Verified landscape (June 2026)

The phenomenon class was opened in the last few months by **De Huang's group**:

- **arXiv:2603.25104** (Huang, Tong, Wang — Mar 2026), *Self-similar finite-time
  blowups with singular profiles of the generalized CLM model*. For the gCLM
  family `ω_t + a u ω_x = u_x ω`, `u_x = H[ω]`: for `a>0`, one-scale self-similar
  blowups with **regular** profiles (new); for `a≤0`, a **two-scale** scenario
  whose **outer profile converges to a singular (unbounded) function**.
  Explicit singular profiles constructed for `a<0`; the `a=0` (CLM) outer-profile
  convergence is **proven** rigorously; inner traveling wave existence via a
  fixed-point method.
- **arXiv:2604.01868** (Chen, Huang, Li — Apr 2026), *Novel self-similar
  finite-time blowups with singular profiles of the 1D Hou–Luo model and the 2D
  Boussinesq equations: a numerical investigation*. A **two-stage** blowup
  (local `L∞` blowup at `T̃`, then weak continuation to a local `Lᵖ` blowup at
  `T>T̃`), with an **unbounded `(X−1)^{−1/2}` self-similar profile** and measured
  exponents `(c_l=2, c_ω=−1, γ=2)`, `1/√X` far-field decay. **Numerical only —
  no rigorous existence of these singular profiles.**

**The gap.** Existing CAP frameworks (Chen–Hou, *Stable nearly self-similar
blowup … II: rigorous numerics*) certify **smooth** profiles. The new **singular
(unbounded)** profiles have, for the Hou–Luo / 2D-Boussinesq case, *no* rigorous
existence proof and (for the whole family) *no published linear-stability
spectrum*. Stability governs which profile is dynamically realized and how hard
a CAP is (instability order controls viscous-perturbation feasibility — cf.
DeepMind arXiv:2509.14185).

**Honest scoop risk.** De Huang's group is driving this hard and owns the CAP
machinery. The robust, hard-to-scoop contributions are (1) **independent
numerical confirmation** of a ~2-month-old single-group claim, and (2) the
**linear-stability spectra**, which the construction papers do not report.

## Candidate targets (the cluster — one solver serves all)

| | target | open? | scoop risk | CAP cleanliness |
|---|---|---|---|---|
| **A** | Hou–Luo singular `(X−1)^{−1/2}` profile: confirm + stability spectrum + desingularized profile-ODE → CAP existence | profile numerical-only; stability unpublished | **high** (freshest, Huang's turf) | needs singularity-adapted basis |
| **B** | gCLM `a`-family: singular-profile transition + **stability spectra λ(a)** across `a` | profiles known (a<0 explicit, a=0 proven); **stability open** | medium | clean 1D profile ODE |
| **C** | De Gregorio `a=1`: spectral gap of the linearization about `A·sin` (regularity-window sub-target) | open; headline near equivalent-to-conclusion | low | cleanest CAP target |

The headline-regularity slices (De Gregorio "all smooth data", supercritical SQG)
are **off-limits for solo closure** — stay on profile / stability / CAP
sub-targets and hand the closure leap off.

## Validation gates (no result trusted until these pass)

1. **CLM exact-solution gate — `validate/validate_clm.py` — PASSED.**
   The engine (`a=0, ν=0`, `ω₀=sin x`) reproduces the closed form
   `ω = 4ω₀/((2−t·Hω₀)² + t²ω₀²)` to **2×10⁻¹⁰** relative `L∞` error, with
   measured `t*=2.000000`, `x*=π` exact. Spatial error saturated below the RK4
   time-floor at N=1024/2048/4096.
2. **Hou–Luo physical-DNS — PARTIAL (c_ω confirmed).** `src/houluo.py` reproduces
   finite-time focusing blowup with **`c_ω=−1` to <0.4%** (matches 2308.01528;
   see `docs/findings_houluo_physical.md`). `c_l` is **not** obtainable from
   physical DNS — the self-similar core `(T−t)^{c_l}`, `c_l>1`, goes sub-grid
   almost immediately. ⇒ need dynamic rescaling.
3. **Singular steady state (Thm 2.3) — CERTIFIED.** `validate/verify_singular_steady.py`:
   closed-form `U(X)=2√(1−X)−2 (X<1), −2 (X>1)` ⇒ **`U(1)=−2` exactly**, steady
   residual of eq 2.4 = **2e-16** on X>1; ℝ-Hilbert tool validated vs Lorentzian
   to 1.4e-6. Independent confirmation of arXiv:2604.01868 Thm 2.3. The ~15%
   error in the *numeric* Hilbert of the singular profile is far-field truncation
   of the `1/√X` tail → mapped grid / desingularization needed for stability.
   See `docs/findings_singular_steady_verified.md`.
4. **(next) Linear-stability spectrum** of `(Ω̄,Θ̄)` — Conjecture 2.4, the open
   prize. Linearize eq 2.4 (closed-form `Ū` removes one singular ingredient) on a
   mapped infinite grid; compute the spectrum / unstable-mode count.
5. **(parallel) Dynamic-rescaling DNS** from degenerate ICs → two-stage
   convergence to the singular profile, **HARD STOP** if exponents drift under
   refinement. Smooth profiles `c_l∈(2,4.53)` for context.

## Layout

```
src/gclm.py          unified gCLM/De Gregorio engine (Hilbert, 2/3-dealias, RK4, adaptive dt)
src/clm_exact.py     CLM closed-form solution + blowup-time locator (self-contained derivation)
validate/            falsification gates (gate #1 = CLM exact, PASSED)
results/  docs/       run outputs and notes
```

## Method

Measure → mechanism → remeasure. Dynamic-rescaling measures the singular-point
exponents; a desingularizing change of variables turns the singular profile into
a bounded fixed point and states the exact functional-analytic obstruction the
smooth-profile CAP framework fails at; interval-arithmetic / Lean closes it or
hands off the obstruction at maximum resolution.
