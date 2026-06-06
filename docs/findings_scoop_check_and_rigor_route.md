# Finding: scoop-check (our stability contribution is NOT scooped) + the rigor route is
# ANALYTICAL weighted-energy finite-codimension stability, not a spectral CAP

Date 2026-06-06. A multi-angle adversarial literature sweep (workflow: 5 search angles, 8
candidates, per-candidate skeptical verification against arXiv abstracts) to answer: is the
rigorous existence and/or stability of the singular Hou–Luo / gCLM / 2D-Boussinesq
self-similar profiles already published?

## Verdict

- **STABILITY (our contribution): OPEN.** No paper provides ANY stability analysis (linear,
  nonlinear, or spectral) for these singular (unbounded) self-similar profiles. Our numerical
  stability spectra (the whole `validate/` program) is un-scooped and is the cleanest,
  highest-value contribution.
- **RIGOROUS EXISTENCE (CAP): open / narrow.** The only rigorous existence touching these
  profiles is INSIDE the base paper `arXiv:2603.25104` and is narrow: (i) `a=0` outer-profile
  convergence to an explicit singular function, (ii) an inner traveling-wave fixed-point. A
  general-profile CAP of existence is NOT published — but it is the **scoop-exposed** part
  (De Huang's group owns the numerics + the special-case rigor).

## Competitive landscape (verified, all fetched from arXiv)

| paper | who / when | what it actually does | relation |
|---|---|---|---|
| 2604.01868 | Chen, De Huang, Li — Apr 2026 | Hou–Luo + 2D Boussinesq singular profiles, NUMERICAL ONLY | base (target object) |
| 2603.25104 | De Huang, Tong, Wang — Mar 2026 | gCLM singular profiles; rigorous only for `a=0` outer + inner TW fixed-point | base; narrow rigor |
| 2605.15130 | Jiajie Chen — May 2026 | RIGOROUS existence + **finite-codimension stability** of 3D Euler `C^{1,1/3−}` profiles | feasibility template; wrong model + smooth/Hölder profile |
| 2605.15149 | Jiajie Chen — May 2026 | rigorous `C^∞` 1D limiting profiles for 3D Euler | feasibility; wrong model/profile |
| 2603.11301 | Hou, Qin, Sire, Wu — Mar 2026 | rigorous self-similar blow-up for 1D reduction of generalized SQG | different model; no stability |
| 2509.14185 | Wang, Lai et al. — Sep 2025 | ML/numerical unstable singularities (3D Euler-bdry, IPM) | prior art; numerical; wrong models |

None scoop the singular-profile stability. The Jiajie-Chen pair confirms rigorous
existence+stability of self-similar profiles is METHODOLOGICALLY achievable in adjacent
problems.

## The rigor route (reframed by the template)

`arXiv:2605.15130`'s stability is built from **"anisotropic weighted estimates" in a
"low-regularity setting"**, a **"finite codimension stability argument"**, and
**"integration-by-parts along trajectories that exploits the equation twice"** (133 pp). So
the rigorous-stability route in this lineage is **ANALYTICAL weighted-energy + finite
codimension**, NOT a spectral interval-arithmetic CAP.

This reconciles the CAP-scaffold result (`findings_cap_scaffold.md`): the brute-force
symmetric-part / dissipativity energy estimate FAILS for the full non-normal operator (the
nonlocal × convective-weight non-commutation blows up the numerical range). The rigor works
because:
- **finite codimension** = project out the neutral symmetry subspace `{ξ@1, µ, ν}` (exactly
  our weight-excluded directions) — stability is claimed only on the complement;
- **anisotropic weights** = a sharper, possibly direction/region-dependent weight (not a
  single `e^{-2γz}`) tuned so the nonlocal term is dominated;
- **IBP-along-trajectories-twice** = a structural cancellation for the transport+nonlocal
  coupling, the analog of our analytic transport estimate `−2γ‖g‖²_w` but carried to the full
  nonlinear, nonlocal estimate.

The numerics we have (no `Re>0` mode, convective gap, no transient growth, finite-amplitude +
co-evolving relaxation, stress-tested) are the EVIDENCE that such an estimate exists and the
finite-codimension subspace is `{ξ, µ, ν}`; the analytical estimate itself is the leap.

## Status / recommendation

- The un-scooped, COMPLETE contribution is the **numerical stability spectra** — ready to
  consolidate into a writeup (the paper the base works lack).
- The **rigorous stability proof** is a major analytical undertaking (the neighbor is 133 pp):
  the human-leap. Method identified (analytical weighted-energy finite-codimension à la Chen
  2605.15130), obstruction mapped (the sharp anisotropic-weight estimate dominating the
  nonlocal term modulo `{ξ,µ,ν}`), feasibility established (Chen does it for neighbors). This
  is Bryan-led; not a spectral CAP, not interval arithmetic.
- The **existence CAP** is open but scoop-exposed — lower priority than the stability writeup.

Recommended: (1) consolidate the numerical stability-spectra contribution for write-up
(un-scooped, done); (2) pursue the analytical finite-codimension stability proof as the
Bryan-led rigor project, using Chen 2605.15130 as the template and our numerics as the
quantitative guide.
