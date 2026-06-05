# HANDOFF — singular-blowup-1d

Read this + `README.md` + `docs/findings_*.md` first. The next step is named at
the bottom ("THE NEXT CHIP"). Repo: github.com/Brsanch/singular-blowup-1d (private).
Started 2026-06-05. Single-researcher + AI; tight measure→mechanism→remeasure loop.

## The bet (one paragraph)

Certify *singular* self-similar finite-time blowup of 1D nonlocal active-scalar
models (CLM / De Gregorio / gCLM / **Hou–Luo**). This INVERTS the NS/SQG
regularity wall: instead of ruling out blowup on the true equation (→ Type-I /
KNSS-Liouville rigidity, equivalent-to-conclusion), we *certify* blowup on toy
models where blowup is the answer. 1D ⇒ extreme resolution cheap; the profile
ODEs are the cleanest CAP/Lean targets. Chosen over SQG-patch / OSW / gCLM-
classification by a 27-candidate adversarial scout (`docs/candidate-scan.md`).

## Verified landscape (the two anchor papers, both real, checked)

- **arXiv:2603.25104** (Huang, Tong, Wang, Mar 2026): gCLM family `ω_t+a u ω_x=u_x ω`,
  singular profiles for `a<0` (explicit), `a=0` (CLM) rigorous, `a>0` regular.
- **arXiv:2604.01868** (Chen, **De Huang**, Li, Apr 2026): 1D Hou–Luo + 2D
  Boussinesq, two-stage `L∞→Lᵖ` blowup, **unbounded `(X−1)^{−1/2}` profile**,
  exponents `c_l=2, c_ω=−1, γ=2`. NUMERICAL ONLY — no rigorous existence, **no
  published stability spectra**. ← our hard-to-scoop target.

**Scoop risk:** De Huang's group owns the CAP machinery and is driving this hard.
Robust contributions: independent confirmation (done) + stability spectra (open).

**Target (Bryan picked):** "confirm + stability spectra, Hou–Luo lead."

## The Hou–Luo model + the certified base state

Model (`x∈ℝ`): `ω_t + u ω_x = θ_x`, `θ_t + u θ_x = 0`, `u_x = H[ω]`.
Dynamic-rescaling system eq.(2.4) (steady state = the profile; linearization =
the stability operator):
```
Ω_τ + (U + c_l X)Ω_X = c_ω Ω + Θ_X
Θ_τ + (U + c_l X)Θ_X = (c_l + 2c_ω)Θ,   U_X = H(Ω), U(0)=0
```
**Certified singular steady state (Thm 2.3), closed form:**
```
Ω̄(X) = 1{X>1}/√(X−1),  Θ̄(X) = (π/2) 1{X>1},  c_l=2, c_ω=−1,
Ū(X)  = 2√(1−X)−2  (X<1),  −2  (X>1)   ⇒  U(1)=−2 exactly,
ā(X) := Ū+2X = 2(X−1) on X>1.
```

## State of play (all committed + pushed)

| result | status | where |
|---|---|---|
| gCLM engine vs CLM exact | ✅ 2e-10 | `validate/validate_clm.py` PASS |
| Hou–Luo physical DNS `c_ω=−1` | ✅ <0.4% | `findings_houluo_physical.md` |
| Thm 2.3 singular profile certified (`U(1)=−2`, residual 2e-16) | ✅ | `validate/verify_singular_steady.py` PASS, `findings_singular_steady_verified.md` |
| dynamic-rescaling solver + "profiles are saddles" | ✅ | `findings_dynamic_rescaling.md` |
| stability spectrum (frozen `μ=ν=0`): resolved sector stable; leading mode = free-boundary artifact | ✅ evidence-for-stability | `findings_stability_spectrum.md` |

## File map

```
src/gclm.py            gCLM/De Gregorio engine (CLM-exact-validated)
src/clm_exact.py       CLM closed form (self-contained Riccati derivation)
src/houluo.py          physical-space Hou-Luo solver (c_omega=-1)
src/rescaled.py        R-line Hilbert (validated 1.4e-6) + singular profile + steady residual
src/dynamic_rescaling.py  eq (2.4) time solver (saddle finding)
src/stability.py       2N x 2N linearized operator (frozen mu=nu=0), closed-form U_bar
validate/              all gates + probes + the spectrum convergence study
docs/                  candidate-scan + 4 findings docs
```
Run any gate: `python3 validate/<name>.py`. numpy/scipy present. Pure-Python, 1D,
cheap (no GPU). No kernel-panic risk (no du/find/lake here).

## Dead ends (do NOT redo)

1. Forward time-marching of eq (2.4) to reach SMOOTH profiles — they are saddles
   (nonlinear via `U=H(Ω)`); constraint-`c_ω` runs away (+39), fixed `c_ω=−1`
   decays to 0. Use Newton/BVP if smooth profiles are needed (they are NOT needed
   for the singular-profile stability deliverable — singular profile is closed-form).
2. NAIVE regularized eigenvalues of the singular operator — leading modes are
   δ-divergent grid-scale EDGE artifacts (`−0.13 → +24+36i` as δ:0.14→0.035,
   `x*≈1`, hi-k~0.45). Do NOT report these as instabilities.
3. `c_l` from physical DNS — the self-similar core `(T−t)^{c_l}`, `c_l>1`, goes
   sub-grid immediately. Use dynamic-rescaling coords.

## THE NEXT CHIP — desingularized free-boundary stability eigenvalue problem

The definitive Conjecture-2.4 verdict needs the operator handled at the support
edge `X=1`. The frozen frozen-coefficient resolved spectrum is already stable;
the open piece is the edge mode. Formulation (derived, ready to implement):

- **Edge-fixed frame** `Y = X − εξ(τ)`, `ξ` = edge displacement (scalar DOF).
  X-frame perturbations: `ω = G − ξ Ω̄'`, `θ = H − ξ Θ̄'`, with `G,H` REGULAR
  shape perturbations. The `−ξΘ̄' = −ξ(π/2)δ(Y−1)` reproduces the delta-source
  (edge motion), removing it from the shape equation.
- **Velocity is well-defined:** `u = H(ω) = H(G) − ξ H(Ω̄')`, and
  `H(Ω̄') = (H Ω̄)' = Ū_XX` (H commutes with ∂_X), which is `0` for `X>1` and
  `−½(1−X)^{−3/2}` for `X<1` — finite as a function on the support side.
- **Edge compatibility:** the `(X−1)^{−3/2}` coupling `u·Ω̄_X` is finite only if
  the relevant edge combination (≈ `u(1)+μ` paired with `ξ`) cancels; this yields
  the `ξ`-ODE. Work in the weighted space where perturbations scale like the
  profile `(X−1)^{−1/2}` (NOT its derivative `(X−1)^{−3/2}`); the edge condition
  is effectively `Hω(1)→0`.
- Build `eig` of the `(G, H, ξ)` operator (closed-form `Ū` removes one singular
  ingredient for free). **Then add the `μ,ν` scaling/symmetry modes** (the
  normalization directions) and project them out before reading the asymptotic-
  stability verdict.
- **Validation targets:** (i) a known scaling/symmetry mode should appear at its
  predicted eigenvalue; (ii) eigenvalues converge under mesh + weight refinement
  AND are edge-resolved (low hi-k); (iii) cross-check against a complementary
  NONLINEAR test — dynamic-rescaling DNS (`src/dynamic_rescaling.py`) from a
  DEGENERATE IC (`ω_x(0)=θ_xx(0)=0`) with a stabilized normalization, watching
  two-stage convergence to `Ω̄` (the paper's own method). HARD STOP if exponents
  drift under refinement.

⚠️ This is the one place a careless discretization yields a plausible-but-WRONG
spectrum. Validate every step; do not report an unstable λ whose eigenmode is
edge-spiked/grid-scale or δ-dependent.

## Bryan's working preferences for this repo
Single-agent serial; local-verify-primary; commit each gated milestone; never
report a spectrum/instability without a convergence + resolution check
(no-overselling); state the obstruction at maximum resolution and hand off the
non-interpolative analytical leap rather than grinding past it.
