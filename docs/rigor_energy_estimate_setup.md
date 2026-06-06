# Rigor setup: the finite-codimension weighted-energy stability estimate (chip #3 entry point)

Date 2026-06-06. The maximum-resolution analytical setup for a RIGOROUS stability proof of
the singular Hou–Luo self-similar profile (Conjecture 2.4), following the route identified
in `findings_scoop_check_and_rigor_route.md` (analytical weighted-energy + finite codimension,
à la Jiajie Chen arXiv:2605.15130 — NOT a spectral CAP). This document states the exact
inequalities to prove. It does not prove them; the sharp coercivity estimate is the leap.

## 1. Operator and variables (established)

Moving (edge-fixed) frame, frozen `μ=ν=0` core (`findings_edge_condition_derivation.md` §4),
base `Ω̄=(Y-1)^{-1/2}`, `Θ̄=(π/2)1{Y>1}`, `c_l=2`, `c_ω=-1`, `ā=Ū+c_l Y=2(Y-1)` on `Y>1`.
Perturbation `Ω=Ω̄+g`, `Θ=Θ̄+h`. Bounded edge-resolved variables: `z=ln(Y-1)`, `G=(Y-1)^{1/2}g
= e^{z/2}g`, `∂_Y = e^{-z}∂_z`. The linearized shape system is, in `(G,h)`:

  `G_τ = M_GG G + M_Gh h`,  `h_τ = M_hh h`,  with the scalar edge ODE `ξ̇=(c_l-1)ξ+V_g(1)+μ`.

where (verified; `src/stability_edge_z.py`, and the transport+`c_ω` algebra collapses cleanly):

  `M_GG = -2 ∂_z + F`,   `F = ½ Q_op`,   `Q_op G = δU/(Y-1)`,  `δU(Y)=∫_1^Y H(g)\,dY'`,
  `M_hh = -2 ∂_z`,       `M_Gh = e^{-z/2} ∂_z`   (the `+h_Y=Θ_X` coupling).

Note the transport and the `c_ω=-1` source COMBINE to the clean `-2∂_z` (the `+G` from
`-ā g_Y` cancels the `-G` from `c_ω g`). `F` is the nonlocal forcing `-[∫_1^Y H(g)]Ω̄_Y`
written in `G`; it is REGULAR at the edge (`~ ½ H(g)(1)(Y-1)^{-1/2}`, `findings_edge_condition`
§4), i.e. `Q_op G → ½ H(g)(1)` finite as `z→-∞`.

## 2. Weighted energy framework

Anisotropic weight `W(z)=e^{-φ(z)}`, `φ'>0` (decay toward the edge `z→-∞`; outflow at
`z→+∞`). Inner product `⟨u,v⟩_W=∫ W u v\,dz`, norm `‖u‖_W²=∫ W u²\,dz`. Admissible `G,h`
vanish at both ends of the (mapped infinite) domain: `G(-∞)=0` is the **inflow/Robin edge
condition** (drop the homogeneous `(Y-1)^{-(λ+1)/2}` branch, `findings_edge_condition` §5),
and decay at `+∞`. Energy identity:

  `½ d/dτ (‖G‖_W² + ‖h‖_W²) = ⟨G, M_GG G⟩_W + ⟨G, M_Gh h⟩_W + ⟨h, M_hh h⟩_W.`

## 3. Transport term — CLOSED (exact, analytic)

For any `W=e^{-φ}` and admissible `G`:

  `⟨G, -2∂_z G⟩_W = -∫ φ'(z) e^{-φ} G² dz  -  [e^{-φ}G²]_{-∞}^{+∞}`
                 `= -∫ φ'(z) W G² dz  -  (W G²)|_{+∞}`        (inflow term `(WG²)|_{-∞}=0`)
                 `≤ -∫ φ'(z) W G² dz  ≤  -(inf φ') ‖G‖_W².`

So **the transport is dissipative with local rate `φ'(z)`**; the outflow boundary helps
(sign `≤0`); the inflow boundary vanishes by the edge condition. Same for `M_hh`:
`⟨h,-2∂_z h⟩_W ≤ -(inf φ') ‖h‖_W²`. For `φ=2γz` this is the `-2γ` convective gap
(`findings_spectral_gap.md`). **This term needs no further work.**

## 4. The nonlocal forcing — THE CORE INEQUALITY TO PROVE (the leap)

The only Ω-shape term not yet controlled:

  `⟨G, F G⟩_W = ½ ∫ W(z) G(z) (Q_op G)(z) dz`,  `(Q_op G)(z) = e^{-z} ∫_1^{Y(z)} H(e^{-·/2}G) dY'`.

**Target inequality (Ω-shape coercivity):** find `φ` (with `φ'~2γ` in the bulk, tunable near
the edge) such that, for all admissible `G` in the finite-codimension space of §6,

  `⟨G, F G⟩_W  ≤  ∫ φ'(z) W G² dz  -  δ ‖G‖_W²`   for some `δ>0`.            (★)

Combined with §3 this gives `⟨G, M_GG G⟩_W ≤ -δ‖G‖_W²` (Ω-shape strictly dissipative).

**Why a single exponential weight FAILS (mapped, `findings_cap_scaffold.md`):** with
`φ=2γz`, the weighted numerical range of `F` DIVERGES with edge depth — the similarity
`e^{γ(z_j-z_i)}` makes the nonlocal edge↔bulk coupling unbounded; `C_force=10→29→83`. So (★)
is FALSE for `φ=2γz` as a global field-of-values bound. The proof must instead use:
- an **anisotropic `φ`** with `φ'(z)` enlarged where the nonlocal kernel is strong (near the
  edge), so the local dissipation `φ'(z)W` dominates the local forcing density; and/or
- the **structure of `Q_op`** (it is an averaging operator: `δU/(Y-1)=∫_1^Y H(g)/(Y-1)`,
  bounded by `½H(g)(1)` at the edge), not a generic bounded operator — the field-of-values
  bound is pessimistic; a direct kernel estimate `|⟨G,FG⟩_W| ≤ ∫ K(z,z')|G(z)||G(z')|` with
  the explicit Hilbert/averaging kernel `K` is the honest route.

This is the **non-interpolative leap**: the sharp anisotropic-weight + kernel estimate that
makes (★) hold. Our numerics certify it is TRUE (no `Re>0` mode; convective decay; finite-
amplitude + co-evolving relaxation; stress-tested) — they give the value of `δ` and the
shape of the admissible `G`; the analytic proof of (★) is what is missing.

**Numerical attack on (★) — it CLOSES IN THE BULK (`validate/energy_estimate_search.py`):**
Computing the energy abscissa `λ_max(−diag(φ') + sym_W(½Q_op))` (transport handled exactly
as the diagonal `−φ'(z)`, forcing via its weighted symmetric part) with a **bounded-at-edge
anisotropic weight** (`φ'(z)=2γ·½(1+tanh((z−z₀)/s))`, ramping `0` at the edge → `2γ` in the
bulk; `γ=0.4, z₀=−6, s=2`):
- The global abscissa is `≈+40`, BUT the worst mode sits AT the edge node (`z*≈−11.98`) and
  GROWS with `N` (`35.7→38.2→40.2→41.9` for `N=700→2000`) — the discretized `(Y-1)^{-1/2}`
  edge singularity, not a converged field-of-values.
- **Excluding a thin edge layer `z<z_cut`** (where `G→0` by the Robin condition) the BULK
  abscissa is **negative**: `−0.019` (`z_cut=−10`), `−0.106` (`−8`), `−0.41` (`−6`), `−0.65`
  (`−4`) — increasingly dissipative inward. So the nonlocal forcing IS dominated by the
  anisotropic transport dissipation in the bulk, i.e. **(★) holds away from the edge.**

So (★) reduces to a **matched edge-layer estimate**: a thin layer near `z=−∞` (`Y→1`) where
`G→0` and the forcing is regular (`~½H(g)(1)(Y-1)^{-1/2}`), to be closed by explicit edge
asymptotics + the Robin condition (Chen's region-dependent / "anisotropic" edge treatment),
matched to the bulk coercivity above. That edge-layer estimate + the §5 `H¹` coupling are the
**reduced** leap.

## 5. The coupling term — needs a weighted HIGHER-ORDER energy (structural)

`⟨G, M_Gh h⟩_W = ∫ W G e^{-z/2} ∂_z h\,dz`. By Young, `≤ ε‖G‖_W² + (4ε)^{-1}∫ W e^{-z}(∂_z h)²dz`.
The `∂_z h` term is NOT controlled by the `L²` `h`-energy (`M_hh` only bounds `‖h‖_W`). Hence
the `L²` energy does NOT close: one needs a **weighted higher-order (e.g. `H¹_W`) energy** for
`h` (and `G`), differentiating the equations and re-running §3–§4 on `∂_z G, ∂_z h`. This is
exactly Chen's **"anisotropic weighted estimates in a low-regularity setting"** + the
**"integration-by-parts along trajectories that exploits the equation twice"** (the second
use controls the derivative level). The commutator `[∂_z, F]` and the weight `e^{-z}` factor in
`M_Gh` are the new terms to bound at the derivative level — a second instance of (★).

## 6. Finite-codimension structure (established neutral modes)

Stability is claimed on the complement of the finite-dimensional neutral subspace
`𝒩 = span{ ξ-mode @ λ=1,  μ (dilation),  ν (amplitude) }` (`findings_theta_munu_gap.md`,
`findings_coevolve_*`):
- `ν` (amplitude) generator `↔ G=const` — non-normalizable in `W`, the marginal `k=0` mode;
- `μ`/dilation `↔ X Ω̄_X ~ (Y-1)^{-3/2}` — more singular than the profile, excluded by `W`;
- `ξ` time-translation `↔ λ=c_l-1=1`, the pure-edge mode (confirmed pure-`ξ`, numerically).
These are exactly the three directions the dynamic-rescaling normalization pins
(`c_l=-U(1)`, `c_ω=-c_l/2`, edge-stagnation) and the directions the co-evolving march holds
fixed while the shape relaxes. (★) and §5 are to be proven for `G ⊥_W 𝒩`.

## 7. The estimate to prove (one line)

> There exist an anisotropic weight `W=e^{-φ}` (`φ'>0`, `φ'→2γ` in the bulk) and a constant
> `δ>0` such that for all `(G,h)` in the weighted `H¹_W` space with `G,h ⊥_W 𝒩` and the edge
> condition `G(-∞)=h(-∞)=0`,
>   `⟨(G,h), 𝓛(G,h)⟩_{H¹_W} ≤ -δ ‖(G,h)‖_{H¹_W}²`,   `𝓛=[[M_GG,M_Gh],[0,M_hh]]`.
> Equivalently: prove (★) (and its `∂_z`-level analogue) — the nonlocal forcing `F=½Q_op` and
> the coupling `M_Gh` are dominated by the anisotropic transport dissipation `φ'(z)W`, modulo
> `𝒩`. Then `𝓛` generates a strictly contractive semigroup on `(H¹_W ⊖ 𝒩)` ⇒ asymptotic
> stability of the singular profile modulo symmetries (Conjecture 2.4).

## 8. Status — what is closed, what is the leap

CLOSED (analytic): the operator form (§1); the energy identity (§2); the transport
dissipation `-∫φ'WG²` with both boundary signs (§3); the edge regularity / Robin condition
and the identification of `𝒩` (§5–§6); the regularity of `F` at the edge.
ESTABLISHED (numerical, the quantitative guide): (★) holds with a definite `δ` — no `Re>0`
mode (two discretizations), convective gap `2γ`, no transient growth, finite-amplitude +
fully self-consistent co-evolving relaxation, stress-tested, paper-matched normalization.
REDUCED LEAP (after `energy_estimate_search.py`): the BULK coercivity of (★) is numerically
established with an explicit bounded-at-edge anisotropic weight (bulk energy abscissa `<0`).
What remains:
  (a) a rigorous EDGE-LAYER estimate — a thin layer at `Y→1` where `G→0` (Robin) and the
      forcing is regular; close it by explicit edge asymptotics matched to the bulk coercivity;
  (b) the `H¹_W` coupling closure (§5) — differentiate, control `M_Gh` and `[∂_z,F]` via the
      same bulk+edge split (Chen's IBP-exploiting-the-equation-twice).
These are far smaller than the original global (★). Still analytic (not a spectral CAP, not
interval arithmetic); template Chen 2605.15130; numerics give every constant. HONEST: this is
numerical evidence + a proof STRUCTURE with the hard balance (nonlocal vs dissipation) shown
to close in the bulk — NOT a completed rigorous proof. The edge-layer matched asymptotic is
the genuine remaining analytic work.
