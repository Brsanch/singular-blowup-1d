# desk/ — the convective-weight spectrum: exact structure of `M_GG` in `L²_σ` and the ladder resolvent

Date 2026-07-15. Follows `edge_mellin_symbol.md` (LEM-07/08: point spectrum closed mod `𝒩`;
frontier = essential spectrum / strict gap) and `edge_layer_estimate.md` (energy route dead at
the edge). This doc attacks the LEM-06 essential-spectrum half **in the convective
(edge-unbounded) weight itself** — the weight the prior docs believed was foreclosed by
`cap_scaffold`'s field-of-values divergence — and finds that the *spectral* question there is
not only open but exactly solvable in structure: the operator is unitarily scale-free in the
forcing strength, its energy form diagonalizes in closed form (proving the energy no-go
artifact-free), its point spectrum is empty right of the transport line, and its resolvent has
an explicit 4-Gamma ladder representation, machine-verified end-to-end. The remaining object
for a strict gap is a single named estimate (§8).

Machine gates: `validate/edge_convective_structure.py` (all structural identities below),
`validate/edge_convective_resolvent_scan.py` (the depth/weak-coupling scan, OBS-07).
Everything continuum-exact — **no discrete Hilbert transform anywhere** (the repo's known
edge-bias trap, `edge_layer_estimate.md` §2).

---

## 0. Conventions (pinned once — differ from the ξ-conventions of LEM-07/08!)

- `z = ln(Y−1) ∈ ℝ`; **edge** `z→−∞`, **far field** `z→+∞`.
- Transform `𝒢(s) = ∫ G(z) e^{−sz} dz`, inversion on vertical lines `G = (1/2πi)∫ 𝒢(s)e^{+sz}ds`.
- Weighted space `L²_σ := L²(ℝ, e^{−2σz}dz)` ⟺ `𝒢 ∈ L²` on the line `Re s = σ`
  (two-sided Laplace–Plancherel). `σ > 0` = the **convective / edge-unbounded** family;
  the repo's `L²(e^{−2γz}dz)` is `σ = γ` (spectral-gap runs used `γ = 0.4`).
  Membership heuristics: `L²_σ` demands edge decay `e^{σz}` and *forgives* far-field growth
  up to `e^{σz}`; in particular for `σ > 0` a **far-field constant tail is allowed**, an
  edge constant is not.
- Operator (frozen Ω-sector): `M = −2∂_z + K`, `K = ½Q_op`. LEM-07's ladder in these terms:

  > `K e^{sz} = c(s)·e^{(s−½)z}`, `c(s) = −tan(πs)/(2s+1)` — verified by full continuum
  > quadrature (pv Hilbert in Y-space + causal integral) to 1e−13…1e−15 at `s = 0.15, 0.30, 0.45`.

  Equivalently `K = M_{e^{−z/2}} ∘ T_c` with `T_c` translation-invariant of symbol `c`;
  `K` maps `L²_{σ'} → L²_{σ'−½}` boundedly for `σ' ∈ (−½,½)` (it *lowers* edge decay by ½).
- The homogeneous 4-Gamma solution in z-form (LEM-08's `Ĝ_p` at doubled argument, but note the
  ξ-convention has `s_ξ = −2s_z`-type reflection — do not import LEM-08's pole-to-end
  assignments; re-derived here):

  > **`X(s) = Γ(2s+λ)Γ(2s+2)Γ(s)Γ(1−s)` satisfies `(λ+2s)X(s) = c(s+½)X(s+½)` identically**
  > (machine-verified to 25+ digits over complex `λ, s`).

  `X` is zero-free; line decay `|X(σ+iτ)| ~ e^{−3π|τ|}|τ|^{4σ+Reλ+1}` (Stirling);
  `X(½) = 2πΓ(1+λ)`; `Res₀ X = Γ(λ)`.

## 1. Lemma A (ε-similarity): the convective weight is scale-free in the forcing

On `L²_σ` the translation `(T_h u)(z) = u(z−h)` satisfies `‖T_h u‖_σ = e^{−σh}‖u‖_σ` exactly,
so `U_h := e^{σh}T_h` is **unitary**. Since `−2∂_z` commutes with translations and
`T_h K T_h^{−1} = e^{h/2}K` (from `K = e^{−z/2}T_c`),

> `U_h M U_h^* = −2∂_z + e^{h/2}K` — i.e. **`M` is unitarily equivalent to `−2∂_z + εK`
> for every `ε > 0`** on every pure-exponential weighted space.

Consequences.
1. Spectrum, essential spectrum, resolvent norms, pseudospectra, numerical range, and
   semigroup norms of `M|_{L²_σ}` are all **independent of the forcing strength**.
2. **The `cap_scaffold` paradox dissolves**: on the infinite line the numerical range is
   ε-invariantly what it is (infinite; §3), while any *truncated* computation breaks the
   similarity and manufactures scale-dependent FoV numbers — the measured
   `C_force ∝ e^{−z_min/2}` divergence is the truncation edge re-entering, not a weighted-space
   property that distinguishes weights.
3. **Edge-depth ≡ weak coupling**: for forcing translated to depth `z₀`,
   `‖R_M(λ)T_{z₀}F‖/‖T_{z₀}F‖ = ‖R_{−2∂_z+εK}(λ)F‖/‖F‖` with `ε = e^{z₀/2}`.
   The OBS-07 depth scan *is* the weak-coupling uniformity probe.

## 2. Lemma B (far-field transport line): `{Re λ = −2σ} ⊆ σ_ess(M|_{L²_σ})`

Weyl sequences at `z → +∞`, where the coefficient `e^{−z/2} → 0` kills `K` relative to the
transport (`T_c` is bounded on the line `Re s = σ` since `c` is analytic and `~1/|s|` there;
the kernel tails of `T_c` are `e^{v/2}` left / `(a+bv)e^{−v/2}` right, so mass leakage toward
the edge is exponentially suppressed). The transport symbol on the line is `−2s`, giving the
whole line `{Re λ = −2σ}`. So **the gap can never exceed `2σ`**, and the measured convective
decay rate `2γ` (`spectral_gap.py`) is exactly this far-field line. The natural window is
`σ ∈ (0,½)` (at `σ = ½` the line `Re s = ½` hits the `tan` pole = the `λ=1`/`ξ` symmetry
scale), so the **supremum gap available from this family is `2σ → 1⁻`**.

## 3. Lemma C (energy no-go, artifact-free, with exact mechanism)

Conjugate to flat space: `Ṽ: L²_σ → L²`, `M̃ = −2∂_z − 2σ + K̃`, `K̃ = e^{−z/2}T_{c(σ+·)}`.
Substituting `w = e^{−z/4}u` **diagonalizes the K-quadratic form exactly**:

> `Re⟨K̃u, u⟩_{L²} = (1/2π) ∫ Re c(σ + ¼ + it) · |ŵ(t)|² dt`,  `w = e^{−z/4}u`.

(One line: `⟨e^{−z/2}T u, u⟩ = ⟨(e^{−z/4}Te^{z/4})w, w⟩` and the conjugated symbol is
`c(σ+¼+it)`; legal since `T_c`'s kernel decays `e^{−|v|/2}` > the `e^{|v|/4}` tilt.)

The multiplier profile decides everything. Machine-measured: `Re c(σ+¼) > 0` for **every**
`σ ∈ (0,½)` (values `+7.79, +3.01, +1.40, +0.85, +0.57, +0.43` at `σ+¼ = 0.52…0.76`), with a
**positive band `|t| < t₀(σ)`**, `t₀ ≈ 0.13–0.33`, and `Re c < 0` outside. Hence:

- **No pure-exponential weight is dissipative**: the resonant-exponent tail
  `u = e^{(σ+¼)z}·(window of length L)` concentrates `ŵ` at `t ≈ 0` and gives
  `Re⟨Ku,u⟩_σ/‖u‖²_σ = +0.18, +1.09, +4.04` at `L = 4, 8, 16` (machine gate) — divergent
  positive. Together with `edge_layer_estimate.md` (all bounded-at-edge weights dead), the
  **energy route is now closed for the full exponential family by proof**, replacing the
  truncated `C_force` evidence (which conflated this continuum resonance with the
  discrete-Hilbert edge bias — same `e^{−z_min/2}` rate, different objects).
- Localized bumps are *dissipative* (`Re C(φ) < 0` for all tested modulations, e.g.
  `−0.579` for a unit Gaussian, exact `e^{−a/2}` depth-scaling verified to 8 digits): the
  obstruction is **only** the low-frequency band in `w` — a sharply localized structure a
  modified-energy/hypocoercivity or resolvent argument can aim at.

## 4. Lemma D (the rank-one defect = the firewall as a domain fact)

For strip-regular `G` (i.e. `G ∈ L²_σ ∩ L²_{σ+½}`, transform analytic on `[σ, σ+½]`):

> `K G = (L²_σ-part) − (𝒢(½)/2π)·𝟙`, where `𝒢(½) = ∫G e^{−z/2}dz` and `𝟙 ∉ L²_σ`.

(The `s = 0` simple pole of `c(s+½)`, residue `1/(2π)`, crossed when re-expanding the output
from line `σ−½` to line `σ`; sign machine-pinned: `K(bump) → −𝒢(½)/(2π)` at the edge,
measured `−0.18214` vs predicted `−0.18213`.) Hence on the convective space the generator's
domain carries the constraint

> **`D(M|_{L²_σ}) ⊆ {𝒢(½) = 0}`** — the repo's `B = H(g)(1) = 0` firewall object is the
> *domain condition of the generator*, not an assumption (`½H(g)(1) = −𝒢(½)/(2π)`).

Two exact consequences, both machine-verified:

- **Defect identity.** The inversion `G_p` of `X` on the line solves
  `M G_p = λ G_p − Γ(1+λ)·𝟙` (defect constant `X(½)/2π = Γ(1+λ)`), with far-field constant
  `G_p → Γ(λ)` and edge decay exactly `e^{z}` (measured ratio `e^{−2}` at `z = −30/−28`;
  far value `2.5437` vs `Γ(0.35) = 2.5461` with the predicted `e^{−λz/2}` correction; also
  verified at complex `λ = 0.2+0.9i`). So `X`-inversions are **not** eigenfunctions — the
  naive "the convective weight admits the 4-Gamma modes and the spectrum fills" scenario is
  **false**, killed by the domain defect.
- **Residue⇔value identity.** Any line solution of the resolvent equation satisfies
  `Res₀ 𝒢 = 𝒢(½)/(2πλ)` (from the difference equation at the `s=0` pole; cross-checked
  exactly on `X`: `X(½)/(2πλ) = Γ(1+λ)/λ = Γ(λ) = Res₀X`). Since `Res₀𝒢` = the far-field
  constant of `G`, **the domain condition `𝒢(½) = 0` and "no far-field constant tail" are the
  same condition**: the `𝒩`-marginal structure enters the convective space only through this
  single scalar, and imposing it removes both pathologies at once.

## 5. Proposition E (no point spectrum right of the line)

For `σ ∈ (0,½)` and `Re λ > −2σ`: no `L²_σ` eigenfunction of `M`, **including at the
`𝒩`-values** (`λ = 0`: `𝟙 ∉ L²_σ`; `λ = 1`: the `ξ`-mode has `e^{−z/2}`-type edge growth,
`∉ L²_σ` — the weighted space is the symmetry quotient, as the numerics found).

Proof structure (the LEM-08 rigidity transferred to the line `σ`):
(i) *[edge bootstrap — Lemma E′, the one honest gap]* an `L²_σ` eigenfunction has no
`e^{−λz/2}` homogeneous edge branch (`|e^{−λz/2}|²e^{−2σz}` diverges at the edge exactly for
`Re λ > −2σ` — the threshold *is* the transport line), hence carries the recessive `e^{z}`
edge decay and is strip-regular on `[σ, σ+½]`. This is the continuum form of the repo's
inflow-BC/admissibility step (4-way-derived in `findings_edge_condition_derivation.md`), here
still assumed, not proved — it is the boundary-regularity lemma to write.
(ii) Then `𝒢 = X·P` with `P` ½-periodic; `X` zero-free ⟹ any pole of `P` is a line pole of
`𝒢` ⟹ `P` entire in `w = e^{4πis}` on `ℂ*`.
(iii) Strip control (`H²` of the full period strip = all rays in `w`) + `|X| ~ e^{−3π|τ|}` ⟹
`|P(w)| ≲ |w|^{∓3/4}` near `0/∞` ⟹ Laurent coefficients vanish for `|k| ≥ 1` (`3/4 < 1`) ⟹
`P ≡ const` — LEM-08's integer obstruction, here as a coefficient-decay Liouville.
(iv) Domain (Lemma D): `𝒢(½) = P·X(½) = P·2πΓ(1+λ) = 0` ⟹ `P = 0`. ∎ (mod E′)

## 6. Proposition F (the exact resolvent: ladder + defect closure)

**Closed form.** The variation-of-parameters coefficients collapse to Gamma ratios:

> `Π_n(s) = X(s) / [X(s+n/2)·(λ+2s+n)]` (machine-verified n = 0…11 to 24+ digits), and for
> `F` with entire, exponential-type transform,
> `𝒢₀(s) = Σ_{n≥0} Π_n(s)𝓕(s+n/2)` converges double-factorially and solves the line equation
> `(λ+2s)𝒢(s) − c(s+½)𝒢(s+½) = 𝓕(s)` (machine: rel. residual ~1e−15).

**Defect closure.** `Ĝ := 𝒢₀ − [𝒢₀(½)/X(½)]·X` enforces the domain condition `Ĝ(½) = 0`;
its inversion `G` then satisfies `(λ−M)G = F` **as an operator identity** (machine: PDE
residual at quadrature accuracy ~4e−4 across the line, `λ = 0.3`, `z₀ = −6`), with the
far-field constant of `G` vanishing automatically (§4 residue⇔value identity; measured
`4.5e−5 ≈ 0`). The homogeneous freedom is exactly one constant (`P = const`, by §5(iii));
one freedom ↔ one defect condition: **index-0 closure, no cokernel obstruction on the dense
class**. ⚠️ Subtlety (cost two bugs): `s = ½` is a removable zero×pole point of the ladder —
odd-`n` terms vanish (`c(1) = 0` ⟺ the `Γ(1−s)` pole of `X` at integers) but even-`n` terms
have finite `c(1)c(3/2)`-type limits (`−1/12` at `n=2`), so the naive evaluation at exactly
`½` is wrong (`14.25` vs true `4.30` in the test case) — **and a finite `±δ` limit is *also*
wrong at depth**: the ladder's curvature near `½` grows like the peak `e^{(3/2)n_*\ln n_*}`
(`n_* ≈ e^{|z₀|/3}`), so the `δ²`-error swamps the true value for `z₀ ≲ −12` while being
**dps-stable** (it is an analytic error — it passed a 90/160/260-digit consistency test and
produced spurious `1e12–1e15` gains; the line-equation self-check is blind to the fix constant
since both `𝒢₀` and `X` solve the line equation). The only safe evaluation is the **exact
closed form** `Π_n(½) = X(½)/[X(½+n/2)(λ+1+n)]` summed over even `n`.

So: **`λ ∈ ρ(M|_{L²_σ})` for `Re λ > −2σ` ⟺ the map `F ↦ G` above is `L²_σ`-bounded on the
dense class** (surjectivity then follows by closedness, injectivity by §5). The whole strict-gap
question is now one uniform estimate.

## 7. OBS-07 — the depth / weak-coupling scan (numerics; the falsification probe)

`gain(λ, z₀) = ‖G‖_σ/‖F‖_σ` for `F = e^{−cosh(z−z₀)}`, `σ = 0.4`, λ across
`{0.3, 0.05, −0.2, −0.5, 0.3+1.5i, −0.7}`, `z₀ = −2 … −14`. By Lemma A this scans the
weak-coupling family `ε = e^{z₀/2} → 0`. Cancellation inside the ladder reaches
`~e^{2e^{|z₀|/4}}` (measured 1.2e2 at `z₀=−6`), so the scan runs at dps 90 — the terms are
huge and alternating; **this cancellation is also why any termwise/triangle-inequality bound
must fail at depth, and why the analytic route needs the §8 resummation.**
Gate: `validate/edge_convective_resolvent_scan.py` (quick mode default; `--deep` = this grid).

RESULTS (2026-07-15, dps 90, cancellations up to `10^67` handled with ~23-digit margin;
`gain(0.3, −6) = 0.08173` independently reproduced by the float64 numpy implementation and
PDE-verified there):

| `z₀`   | −2          | −6            | −10     | −12     | −14     |
|--------|-------------|---------------|---------|---------|---------|
| gain   | 0.303–0.370 | 0.0817–0.0852 | 0.01337–0.01350 | 0.00503–0.00505 | 0.001868–0.001871 |

across all six λ (real from `0.3` down to `−0.7` = dist 0.1 from the line, plus `0.3+1.5i`).
Successive-depth log-slopes: `−0.327, −0.454, −0.490, −0.494` → the asymptotic law is

> **`gain(λ, z₀) ∝ e^{z₀/2} = ε`, λ-uniform** — the response to edge forcing *vanishes
> linearly in the weak coupling*; it does not merely stay bounded.

Reading: (i) the resolvent-bound question is decided *favorably on the edge-forcing family
with large margin* — the sup of `‖R(λ)F‖/‖F‖` lives at bulk-scale forcing, where the termwise
bound already gives `≈1.15/dist` (§8); together this is strong (numerical, not proof) support
for `σ(M|_{L²_σ}) = {Re λ = −2σ}` with `‖R(λ)‖ ~ C/dist`. (ii) The mechanism matches Lemma C:
the K-form on edge-localized structures is dissipative with the `e^{|a|/2}`-large coefficient —
the operator *annihilates* edge input. (iii) A genuine singular-perturbation discontinuity:
`lim_{ε→0}‖R_{A+εK}F₀‖ = 0 ≠ ‖R_A F₀‖` — the naive weak-coupling expansion around pure
transport is invalid on fixed data (at depth the ladder is dominated by high `X`-shaped rungs
and `C₀(line) ≈ C₀(½)` to 15+ digits; the transport component is exponentially subdominant).
(iv) ⚠️ Numerics lesson recorded in §6: the first run of this scan produced dps-STABLE garbage
(`1e12–1e15` gains at `z₀ ≤ −12`) from a `±δ` defect-constant evaluation — an *analytic*
error no precision test catches; only the closed-form `Π_n(½)` evaluation is safe.

## 8. The sharpened remaining object (the leap, updated)

The strict gap `Re σ(M|_{L²_σ}) ≤ −2σ + δ`-half of LEM-06 is now exactly:

> **Bound the defect-closed ladder resolvent uniformly on `L²_σ`, `Re λ ≥ −2σ + δ`.**
> Its most compact exact form: with `C₀(s) := Σ_n 𝓕(s+n/2)/[X(s+n/2)(λ+2s+n)]` (which
> satisfies the coboundary relation `C₀(s) − C₀(s+½) = 𝓕(s)/[X(s)(λ+2s)]` exactly),
> the solution is
> **`𝒢(s) = X(s)·[C₀(s) − C₀(½)]`** —
> at depth `C₀` is huge but nearly ½-periodic (its defect from periodicity is the single
> `n=0` term), and the resolvent is `X` times the *deviation* `C₀(s) − C₀(½)`; the analytic
> work is bounding that deviation (the "two-stage cancellation" the numerics exhibit).
> Named route (new this session): the closed form makes `𝒢₀/X` a Barnes-type sum
> `Σ_n 𝓕(s+n/2)/[X(s+n/2)(λ+2s+n)]`; rotate it to a contour integral (resummation kills the
> cancellation problem) and Schur-test the resulting explicit Gamma-ratio kernel. ⚠️ The
> vertical-line rotation **diverges** — the integrand grows like `e^{π|t|/4}` along vertical
> ν-lines (measured slopes 0.47/0.64/0.73 → π/4; `1/X`'s `e^{3π|t|/2}` beats the sin-kernel's
> `e^{−π|t|}` plus the transform's `e^{−π|t|/4}`) — so the contour must be **bent** into the
> right half-plane (Stirling real-direction growth of `X` then gives superexponential decay),
> and the kernel must respect the ladder's internal period-2 sign structure
> (`1/Γ(1−s−ν/2)`'s reflection factor). This is a bounded, classical-flavored computation —
> graded Mellin–Barnes asymptotics — not an open-ended search.

Quantitative anchors for that work: (i) the termwise line-sup constants
`M_n(λ) = sup_τ |Π_n(σ+iτ)|` are tame — `Σ_n M_n ≈ 1.15/dist(λ, line)` uniformly (measured
`1.02/0.909`, `3.75/3.33`, `11.5/10.0` at dist `1.1/0.3/0.1`): the K-dressing of the *termwise*
bound is a ~15% overhead on pure transport, which proves the two-space estimate on the
superexponential-edge class outright and suggests the full `L²_σ` bound carries the same
`C/dist` shape. (ii) The weak-coupling limit is NOT the naive transport resolvent: the
ε-scaled homogeneous solution is `X_ε(s) = ε^{−2s}X(s)` (i.e. `X` translated to the forcing
depth), so the defect-fix term never becomes relatively small — the constrained resolvent
converges to a *constrained* transport problem. The scan (§7) shows the fix in fact
*suppresses* the response on edge data.

Also still open, named: **Lemma E′** (edge bootstrap / boundary regularity of eigenfunctions
and of the maximal domain — the continuum version of the inflow BC); **semigroup generation**
on the constrained domain + Gearhart–Prüss to convert resolvent bounds into
`‖e^{tM}‖_σ ≤ Ce^{(−2σ+δ)t}`; then the `M_Gh` coupling (same `e^{−z/2}`-ladder structure) and
the nonlinear/co-moving leg (unchanged).

Classical anchors (context, not off-the-shelf solutions — our edge end has *no* limit operator;
the ε-similarity replaces that role): exponential weights for convective stability go back to
Sattinger (Adv. Math. 22, 1976); the essential-vs-absolute spectrum frame is Sandstede–Scheel
(Physica D 145, 2000) — their Fredholm-border machinery is ODE-based and does not apply verbatim
to the `½`-shift ladder; resolvent-to-semigroup is Gearhart–Prüss (Engel–Nagel, *One-Parameter
Semigroups*, Thm V.1.11); the Barnes/Mellin factorization toolkit is Gohberg–Krein-adjacent but
the graded `½`-shift with the `𝒩`-defect appears genuinely nonstandard (as LEM-07 §8 flagged).

**What this doc changes in the LEM-06 frame:** the essential-spectrum half was "move the
`Re s = 0` line off the axis ⟹ needs the edge-unbounded weight ⟹ re-triggers cap_scaffold" —
a dead end by energy methods. It is now: the edge-unbounded weight's *spectral* problem is
exactly solvable in structure; the energy obstruction is understood (band-limited positive
multiplier, §3) and bypassed (resolvent, §6); the gap question is one explicit uniform bound
(§8) that numerics can probe artifact-free (§7). The `Re s = 0` non-Fredholm line of the
bounded weight is simply the `σ → 0⁺` end of the family — nothing needs to "move"; one works
at `σ ∈ (0,½)` from the start where the transport line sits at `−2σ`.

## 10. EXECUTED 2026-07-16 — the slit-plane Plemelj representation closes the Schur bound (proof skeleton complete)

Gate: `validate/edge_slitplane_representation.py` (3 checks, PASS). This section executes §8's
program one large step: the bent-contour Barnes idea is **superseded by a cleaner classical
reduction** that avoids ν-contours entirely.

**(a) The jump problem.** Under `w = e^{4πis}` the strip `[σ, σ+½]` maps onto `ℂ*` slit along
the ray `arg w = 4πσ` (the two strip boundaries land on the two sides of the slit), and the
coboundary `C(s) − C(s+½) = Φ(s)`, `Φ = 𝓕/[X·(λ+2·)]`, becomes an **additive Plemelj jump
problem** across the ray. The solution is the **split-regularized Cauchy transform** (plain
kernel on the `w'→0` wing, one subtraction on the `∞` wing — Muskhelishvili class-h
bookkeeping forced by the density growth `|Φ| ~ e^{5π|τ'|/2}` from `1/X`), and slit-plane
Laurent rigidity (growth `< |w|`) gives `C₀ = D + const`, so the defect-closed resolvent is
constant-free:

> **`𝒢(s) = X(s)·[D(e^{4πis}) − D(1)]`, `D` = split-regularized Cauchy transform of `Φ`.**
> Machine-certified ≡ the ladder solution at strip-interior points to **rel 1e−22** (jump
> orientation `D₋ − D₊ = Φ` in the τ'-parametrization). Crucially the formula consumes ONLY
> line values of `𝓕` — it defines the solution for **arbitrary `L²_σ` data**; the
> minus-Hardy wall has an explicit classical solution formula.

**(b) The dressed-Schur closure.** The `L²(line) → L²(line)` bound of `𝓕 ↦ 𝒢` reduces to the
two-weight Cauchy kernel `|X(τ)|/|X(τ')|·ker(τ,τ')/|λ+2s'|` with `ker` the split kernel. The
earlier "corner failure" (input `|τ'| ≫ |τ|`, both negative) evaporates: **it used the
unsubtracted kernel on the wrong wing** — the subtraction supplies exactly the missing
`e^{−4π(τ−τ')}` there. All four quadrants then decay at rate `≥ e^{−π|τ−τ'|}` (the `4π − 3π`
margin — the same ¾-ratio as LEM-08's integer obstruction), machine-checked:
`dressed × e^{+π|τ−τ'|}` = 0.117→0.022 (former corner), 3.4e−4→5.4e−21, 1.3e−16→2.8e−33
(others) — bounded and *decaying*, so the polynomial Stirling prefactors have strict
exponential slack. The diagonal band `|τ−τ'| ≤ 1` has weight-ratio `≤ e^{3π+}` and reduces to
local Riesz-projection boundedness + the `½Φ` Plemelj term; the rank-one defect functional
`𝓕 ↦ C₀(½)` must be bounded via the **ladder sum / split measure**, not via bare line
density of `1/(|X|·|λ+2s|)` — that bare density **grows** as `e^{γ|τ|}` with
`γ≈3π` (gate `edge_schur_envelopes.py`, 2026-07-22; desk prose that claimed
`~e^{−3π|τ'|/2}` decay is **false** for `1/|X|` and must not appear in W1).
λ-dependence: `1/dist` from `(λ+2s')^{−1}`, with `δ`-dependent prefactors as
`Re λ ↓ −2σ+δ` (the `Γ(2s+λ)` pole approaching the line).

**(c) OBS-08 — the high-frequency corner probed and benign.** For
`F = e^{−cosh(z−z₀)}e^{iΩz}` (spectrum at `τ ≈ Ω` — the region OBS-07 never probed):
dressing = gain/transport-gain **rises monotonically to 1** with `|Ω|` at every depth
(z₀=−2: 0.457→0.986→0.9996→1.0000 at Ω=0/−5/−15/−40; z₀=−6: 0.123→0.888→0.9968→0.9999;
z₀=−10: 0.020→0.441→0.959→0.9985; ±Ω symmetric). High-frequency edge forcing feels the
weakened coupling `c(σ+iΩ) ~ 1/(2|Ω|)` and responds as pure transport; the low-frequency
suppression crossover matches §3's band structure. Nowhere in the probed `(z₀, Ω, λ)`-region
does the dressing exceed ~1.

**Status after this section (honest).** The uniform-resolvent-bound statement of §8 now has a
**complete proof skeleton**: representation (machine-certified), off-diagonal Schur (exponent
arithmetic + machine-checked envelopes with slack), diagonal band (classical Riesz), defect
functional (explicit kernel), uniqueness (Prop E mod E′), extension to all of `L²_σ`
(closedness of the maximal operator + bound on the dense class). It is NOT yet a proof:
**(W1)** the rigorous write-up (Plemelj boundary-value theory application, Stirling prefactor
bookkeeping, the operator-identity for `D`-defined solutions distributionally) — bounded,
classical-flavored, numerically de-risked at every step; **(W2)** Lemma E′ (edge bootstrap,
for uniqueness); **(W3)** semigroup generation + Gearhart–Prüss to convert `‖R(λ)‖ ≤ C/dist`
into `‖e^{tM}‖_σ ≤ Ce^{(−2σ+δ)t}`. Then the `M_Gh` coupling and the nonlinear/co-moving leg.
The *linear strict-gap wall has no missing idea left* — it is now technical execution.

## 11. Circularity firewall (covers §§1–10)

Nothing here asserts stability from benign numerics: Lemmas A–D and Prop. F are identities
(machine-verified to 13–30 digits); Prop. E is structural (Liouville + Stirling + the domain
condition), conditional only on the *named* Lemma E′; §7 is labeled numerics; §8 is labeled
open. The firewall object `B = H(g)(1)` is *derived* as the generator's domain constraint
(Lemma D), used nowhere as an assumption. The `𝒩`-quotient enters only through the space
choice (§5), matching `theta_munu_gap.py`'s "the weighted space is the symmetry quotient".
§10 adds no stability assumption anywhere: the representation identity is certified against
the independently-verified ladder; the Schur envelopes are computed from the explicit `X` and
kernel (spectrum never enters); OBS-08 is labeled numerics; and the section's status paragraph
states exactly which steps remain unwritten (W1–W3) — "proof skeleton" is not "proof".
