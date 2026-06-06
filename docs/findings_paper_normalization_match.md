# Finding: anchor papers verified real; co-evolving normalization MATCHES arXiv:2604.01868 (eq 3.2)

Date 2026-06-06. Code: `validate/coevolve_paper_match.py` (PASS). Closes the referee-facing
item of THE NEXT CHIP #2 (match the paper's explicit normalization) and verifies the
foundation against the primary sources.

## Anchor papers verified (anti-fabrication on the foundation)

Both fetched from arXiv and confirmed to exist and match the handoff:
- **arXiv:2604.01868** — "Novel Self-similar Finite-time Blowups with Singular Profiles of
  the 1D Hou–Luo Model and the 2D Boussinesq Equations: A Numerical Investigation,"
  Bojin Chen, De Huang, Xiangyuan Li (submitted 2026-04-02). Two-stage `L∞→Lᵖ` blowup,
  unbounded singular self-similar profile. **Numerical only** (no rigorous existence, no
  stability spectra) — our stability-spectrum contribution remains the hard-to-scoop part.
- **arXiv:2603.25104** — "Self-similar finite-time blowups with singular profiles of the
  generalized Constantin-Lax-Majda model: theoretical and numerical investigations,"
  De Huang, Jiajun Tong, Xiuyuan Wang (submitted 2026-03-26). gCLM with parameter `a`:
  `a>0` regular, `a≤0` two-scale singular, `a=0` rigorous. (Matches the gCLM engine
  validated in `validate_clm.py`.)

## The paper's normalization (eq 3.2), quoted from the source

```
U(1,t) + c_l ≡ 0                 (pins the singularity at x=1)
U_x(0,t) ≡ U_x(0,0) = -1         (fixes the velocity gradient at the origin)
```
Singular profile (Thm 2.3): `(Ω̄, Θ̄, c_l, c_ω) = (1{X>1}/√(X-1), (π/2)1{X>1}, 2, -1)`.
(The paper also gives an alternative set 3.3 with `Ω_x(1,t)≡0` for adaptive meshes.)

## Match (validated)

- **Condition 1 is identical.** Our co-evolving normalization uses `c_l = −U(1)`
  (edge-stagnation) — this IS the paper's eq (3.2) line 1 (`U(1)+c_l=0`). It holds in our
  march by construction: `U(1)+c_l = 0` to machine precision at every τ.
- **Profile satisfies condition 2 and the exponents.** On the exact profile:
  `U(1)=−2.000000`, `U_x(0)=H(Ω̄)(0)=−1.000000`, `c_l=2.000000`, `c_ω=−1.000000`.
- **Condition 2 holds asymptotically along our march.** We used the Θ-balance `c_ω=−c_l/2`
  for `c_ω` (it removes the Θ-amplitude h-forcing instability of the perturbation-split
  scheme; the paper's `U_x(0)=−1` would reintroduce it in *our* discretization, though the
  paper's adaptive full-field scheme handles it). Along our relaxing trajectory
  `U_x(0): −1.072 → −1.032 → −1.011 → −1.00002` → **−1**, i.e. the paper's second condition
  is satisfied asymptotically; the two `c_ω` choices agree at and near the profile.

## Verdict

The foundation is real (both anchor papers verified), and our co-evolving stability march
uses the paper's exact `c_l` normalization (eq 3.2 line 1) and converges to its `c_ω`
normalization (`U_x(0)=−1`), with the paper's singular profile and exponents `(2,−1)`. The
referee-facing normalization match for THE NEXT CHIP #2 is complete. (Honest: our `c_ω`
forcing uses the equivalent Θ-balance form for numerical stability; it agrees with the
paper's `U_x(0)=−1` at/near the profile and the trajectory satisfies both conditions
asymptotically. Convective/weighted, numerical, not a proof.)

## Remaining (chip #3, the rigor frontier)

Pin the precise weighted space / max gap; then interval-arithmetic / Lean (CAP) closure of
"no `Re λ>0` mode of the global nonlocal operator in the weighted space modulo symmetries"
(`findings_edge_condition_derivation.md` §7). Scoop-risky (De Huang's group owns the CAP
machinery); this is the human-leap / rigor step, not a numerical one.
