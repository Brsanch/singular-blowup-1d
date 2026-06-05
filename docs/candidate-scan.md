# Open-problem scan (how this target was chosen)

A multi-agent scout ideated 27 open physics/math problems across 5 lenses
(PDE singularity, rigorous turbulence, discrete/computational, quantum/statmech,
formal-verification+ML), adversarially fact-checked each for current (2024–2026)
openness with web search, and scored tool-fit against the toolkit. 5 candidates
were killed as already-resolved. Top finalists by composite score:

| composite | candidate |
|---|---|
| **8.1** | **CAP of the singular self-similar profile of the 1D Hou–Luo / CCF model (two-stage L∞→Lᵖ, unbounded `(X−1)^{−1/2}` profile)** ← chosen cluster |
| 7.3 | SQG patch boundary curvature singularity (Scott–Dritschel log-spiral) — needs a new Lagrangian contour-dynamics solver |
| 7.0 | OSW gCLM blow-up/regularity transition exponent `a*` — tool-fit 10, but headline near equivalent-to-conclusion wall |
| 6.8 | gCLM stable-profile classification `λ_max(a)` across `a` |
| 6.6 | Formalize-and-extend the Chen–Hou Boussinesq/Euler CAP to a Lean-checked certificate |
| 6.3 | Higher-cumulant FCS anomaly of spin transport in the classical Landau–Lifshitz / Heisenberg chain |
| 6.0 | De Gregorio global regularity on the circle (→ spectral-gap sub-target) |
| 5.5 | Smooth interior (no-boundary) axisymmetric Euler self-similar profile |
| 5.4 | Sharp α-threshold for g-SQG patch singularity (half-plane) |
| 4.6 | gSQG smooth (non-patch) self-similar singularity, α-dependence |
| 4.4 | SQG inviscid-limit anomalous-dissipation universality |
| 4.0 | Bressan mixing conjecture sharp constant |

**Killed as resolved/overtaken:** stability index of neural Boussinesq/Euler
profiles (arXiv:2509.14185); C^{1,α}→smooth Elgindi continuation (Hou–Zhang +
Jiajie Chen May-2026); Onsager-supercritical multifractal mixing; XXZ easy-axis
diffusion-constant divergence; conformal-bootstrap interval certificate (3D
Ising done 2026).

**Why the Hou–Luo / gCLM cluster won:** it is the rare target hitting all three
core capabilities on one object (extreme-1D pseudospectral DNS + mechanism/
exponent extraction + Lean/interval CAP), it *inverts* rather than hits the
documented regularity wall, and a real result (independent confirmation of a
~2-month-old claim) is reachable in days while the CAP/stability prize is a
weeks-to-months stretch. The whole cluster (A/B/C in the README) shares one
solver and one measure→mechanism loop, so the foundation is built before the
precise target is committed.
