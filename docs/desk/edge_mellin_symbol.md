# desk/ — the edge Mellin symbol and the Wiener–Hopf reduction of the stability leap

Date 2026-06-14. Follows `edge_layer_estimate.md` (which proved every weighted-energy
method dead at the edge and rescoped the leap to a non-self-adjoint spectral/Fredholm
bound). This doc **executes** that rescope one decisive step: it computes the edge
Mellin symbol of `M_GG` in closed form, reduces the eigenvalue problem to an explicit
**Mellin difference equation**, and pins the entire Conjecture-2.4 edge obligation to a
single **Wiener–Hopf index**. Symbol + reduction are independently verified (4 methods,
8/8 adversarial verdicts non-circular, mpmath to ~1e-10). The index is **named, not yet
evaluated** — that is the leap.

## 1. The half-line Hilbert is a convolution in `z = log(Y−1)`

`Y−1 = e^z`, `Y'−1 = e^{z'}` turns `H[g](Y) = (1/π)pv∫_1^∞ g(Y')/(Y−Y')dY'` into a
translation-invariant convolution:
`  H[g](z) = (1/π) pv ∫ g̃(z')/(e^{z−z'} − 1) dz',   1/(e^u−1) = ½coth(u/2) − ½.`
So `H` is a **Mellin multiplier** on powers `(Y−1)^p`: `H[(Y−1)^p] = 𝓗(p)·(Y−1)^p`.

## 2. The symbol `𝓗(p) = cot(πp)` (verified, classical)

With `t = e^u`, `𝓗(p) = (1/π) pv∫_0^∞ t^{−p−1}/(t−1)dt = cot(πp)` via the classical
`pv∫_0^∞ x^{a−1}/(1−x)dx = π\cot(πa)` (`0<a<1`; convergence strip `−1<p<0`).

- `𝓗(−½) = cot(−π/2) = 0` reproduces the continuum identity `H[(Y−1)^{−1/2}] = 0`
  (Gradshteyn–Ryzhik 3.222.2) — the profile self-annihilation that anchors everything.
- **Independent verification (workflow `wf_0db34a7c`, 4 methods):** keyhole/residue contour;
  a folding `s→1/s` quadrature needing no cot-table (matched to ~1e-12 across the strip);
  the Fourier-axis form `m(iξ) = −i\coth(πξ)` (`Re ≡ 0`, exposing the `−½` constant tail as a
  rank-one zero-mode piece); and the Erdélyi/Tricomi one-sided (airfoil) Hilbert-of-powers
  formula `(1/π)pv∫_0^∞ t^p/(x−t)dt = \cot(πp)x^p`. All agree; the strip `−1<p<0` is a **true
  convergence boundary** (the truncated integral diverges outside it while `cot` stays finite —
  the out-of-strip value is genuine analytic continuation = the Wiener–Hopf content).
  Numerical gate: `validate/edge_mellin_symbol.py`.

## 3. The forcing is a `−½` Mellin ladder (verified)

`Q_op` is a causal convolution (`e^{−u}1{u>0}`, multiplier `1/(1+p)`) composed with the
multiplications `g=e^{−z/2}G`, `dY'=e^{z'}dz'`, `/(Y−1)=e^{−z}`. For `G=(Y−1)^ζ`:
`  (½Q_op)[(Y−1)^ζ] = −[\tan(πζ)/(2ζ+1)]·(Y−1)^{ζ−1/2},   c(ζ):=−\tan(πζ)/(2ζ+1).`
The output exponent is `ζ−½` — the forcing **inherits the profile amplitude `(Y−1)^{−1/2}`**
(this is *exactly* why every energy method dies at the edge, `edge_layer_estimate.md`). So
`M_GG = −2∂_z + ½Q_op` is **NOT diagonalized by single powers**: it is a Mellin **difference**
operator with step `−½`. Verified exactly by sympy composition tracing `build_blocks_z`, and
numerically in-strip (`ζ=0.3`: predicted `−0.860`, measured `−0.839`).

- Integer `ζ`: `\tan(πζ)=0 ⟹ c=0` — **regular** edge modes, `λ=−2ζ ∈ {0,−2,−4,…} ≤ 0`.
- `ζ=½` (`λ=−1`): `\tan(π/2)=∞` — the forced-constant **resonance** = the `+40` field-of-values
  (OBS-03, κ-crushable, on the g-unbounded branch).
- `ζ=−½` (`λ=+1`): an **order-2 pole** `c(ζ) ~ +1/(2π(ζ+½)²)` (both `tan` and `1/(2ζ+1)` blow
  up) — exactly the time-translation symmetry mode in `𝒩`, removed by quotient.

## 4. The eigenvalue problem IS a Mellin difference equation

Expand on the ladder `G ~ Σ_n a_n e^{(s_0 − n/2)z}`. `M_GG G = λG` rung-by-rung gives the
**first-order linear Mellin difference equation of step ½**:
`  (2s + λ)·a(s) = c(s+½)·a(s+½),    c(w) = −\tan(πw)/(2w+1),   s = −λ/2 the leading exponent.`
Its down-ladder Wiener–Hopf solution is `a(s) = [½-periodic prefactor]·∏_{n≥0} c(s+½+n/2)/(2(s+n/2)+λ)`.
`λ` is an `L²_W` point eigenvalue (mod `𝒩`) **iff** `a` admits a Wiener–Hopf factorization
`a = a_+ a_-` holomorphic and `L²` on the critical line, with no pole crossing the strip between
the admissible **edge** indicial line (`z→−∞`, `Re s=−λ/2`) and the admissible **far-field**
indicial line (`z→+∞`). Both ends obey the *same* difference equation; they differ only in the
analyticity strips (`cot`: `−1<p<0`; causal factor `1/(1+s)`: `Re s>−1`) — so this is genuinely
**one** Wiener–Hopf problem on two strips, not a renamed stability hypothesis.

## 5. The remaining object (maximum resolution — the leap)

> **Compute the Wiener–Hopf index / winding number of the factorization of the symbol
> `[(−2∂_z − λ) + ½Q_op]` across the strip matching the edge indicial line to the far-field
> line, and show it is `0` for all `Re λ ≥ 0` modulo `𝒩 = span{ξ@λ=1, μ, ν}`.**
> Equivalently: the forced-constant amplitude `B = a(0) = H(g)(1)/(2λ)` produced by the
> down-ladder product — "no `Re λ≥0` spectrum" ⟺ "the product cannot make `a(0)=B` consistent
> with `L²` on **both** lines for any `Re λ≥0`, mod `𝒩`."

Regime map in the `ζ = −λ/2` plane:
- `Re λ ∈ [0,1)` ⟺ `ζ ∈ (−½,0]` — **inside** the `cot` strip, symbol regular. The down-ladder
  product shows **no automatic pole** forcing an eigenvalue (encouraging) — but the two-ended `L²`
  connection is still required; the product alone does not decide it.
- `λ = 1` (`ζ=−½`) — on the strip boundary AND the order-2 pole = the symmetry mode, removed by `𝒩`.
- `Re λ ≥ 1` (`ζ ≤ −½`) — **outside** the strip; needs the **analytic continuation** of the symbol.
  This is the decisive regime and the genuine Wiener–Hopf content.

## 5.1 The `ξ=(Y−1)^{1/2}` recurrence, and what the eigenvalue count does (and does NOT) show

**Cleaner form (sympy-verified exact, workflow `wf_e6371da0`).** In `ξ=(Y−1)^{1/2}` transport is
the Euler operator `−2∂_z = −ξ∂_ξ` and the `−½` ladder becomes a **unit shift**:
`M_GG = −ξ∂_ξ + ξ^{−1}c(½ξ∂_ξ)`, and `M_GG G=λG` with `G=Σ a_m ξ^m` gives the unit recurrence
`  a_{m+1} = (m+λ)(m+2)\tan(πm/2)·a_m       [= LEM-07's difference eqn at m=2ζ].`
Every step reproduced independently by sympy, zero discrepancy. Since `cot(πp)` is exact on each
power, this recurrence **is** the full operator inside the strip.

**Argument-principle count = 0 — robust IN-STRIP corroboration, NOT the WH index**
(`validate/edge_argprinciple_count.py`): `#eig(M_GG, edge_bc=True)` with `Re λ ≥ {0.05,0.25,0.5,1}`
is `0` across `N∈{1400,2000}`, `z_min∈{−12,−14,−16}`, and under a g-mass penalty; the contour
winding `(1/2πi)∮ tr((M−λ)^{−1})dλ` over `{0.1≤Re≤6,|Im|≤8}` is `0` (`Im ~1e-15`). ⚠️ **This is
corroborating evidence, not the leap** — three caveats (verified):
1. **Different winding.** (C) winds the spectral parameter `λ` and counts discrete matrix
   eigenvalues; the WH index winds the symbol **exponent `s`** across the strip. "(C)=0 ⟹ WH
   index=0" is a **category error**.
2. **Inherits the firewall.** `edge_bc=True` drops the `z_min` node = enforces admissibility
   (`B=0`) by hand, so "0 unstable eigenvalues" *restates* the firewall, not independent evidence.
3. **Cannot reach `Re λ≥1`.** The discrete Hilbert violates `H[(Y−1)^{−1/2}]=0` at the edge node by
   `−0.445·e^{−z_min/2}` (reproduced to 3 sig figs, growing with depth) — the decisive out-of-strip
   region is uncertifiable numerically.

**Termination is suggestive, not a proof.** The ratio grows super-geometrically
(`|a_{k+1}/a_k| ~ 0.7 k²`; `|a_k| = 1, 1e15, 1e67` at `k=0,10,30`), so naive non-terminating series
have zero radius → only terminating chains (`λ=−2ζ≤0`) are formal `L²`. BUT the first-order half-step
equation carries a **½-periodic prefactor freedom** `A→A·P`; the line-growth argument kills
nonconstant `P` only at the **edge** end, and a genuine eigenfunction must decay at BOTH the edge
(`cot` strip `−1<p<0`) and far-field (causal `Re s>−1`) strips — the two-ended matching IS the index.
And "termination ⟹ `λ≤0`" silently uses the firewall (`L²_W ⟹ B=0`) to drop the `λ=+2,+4`
homogeneous branches the indicial equation alone permits.

**Net (status unchanged):** cleanest recurrence form + robust in-strip corroboration + a sharper
statement of the remaining object — but **no advance past `reduced-to-sharp-WH-factorization`**.
(Minor: `\tan(πw) → +i` as `Im w→+∞`, so the growth is super-geometric `~m²`, not literally factorial.)

## 5.2 The Mellin-line structure — and why the winding is a tautology, not a certificate

A direct run at the factorization (continuum mpmath, immune to the discrete-Hilbert bias;
verified workflow `wf_2d8a4427`). It ruled out the lattice shortcut, gave the explicit Γ-form,
and **corrected an overclaim** — but did NOT advance the status.

**Lattice solutions are edge-divergent for ALL λ ⟹ the eigenfunctions are Mellin-line.** The
terminating-up power series `G(ξ)=Σ a_m ξ^m` (top power `M=−λ` from `R(M)=0`, infinite down-tail)
blows up at the edge `ξ→0` like `exp(c/√ξ)` with `c>0` **essentially λ-independent** (checked
`λ=0.7,2.5,0.5+1.5i,−0.7`; the same form for both signs of `Re λ`, so it does NOT discriminate
stable/unstable). So no power-series/lattice solution is `L²_W` off `Re=0` — "termination" is the
wrong object; the genuine eigenfunctions are **Mellin-line** integrals `Ĝ(s)`.

**The Γ-form.** A particular solution of `Ĝ(s+1)=R(s)Ĝ(s)`, `R(s)=(s+λ)(s+2)tan(πs/2)`, is
`Ĝ_p(s)=Γ(s+λ)Γ(s+2)·τ(s)` with the **tan-Barnes** `τ(s+1)/τ(s)=tan(πs/2)` (double-gamma-type).

**The winding — corrected, and it is the indicial-root tautology, NOT the Fredholm index.**
`W(λ)=(1/2π)[arg R(iτ)]_{−∞}^{+∞}` decomposes **exactly** per factor as
`sign(Re λ)/2 + ½(s+2) + ½(tan)`, giving `W=+3/2` (`Re λ>0`, constant) and `+½` (`Re λ<0`) — a jump
of **exactly +1** at `Re λ=0` (one `s=−λ` zero crossing the line). *(An earlier `−0.5`/jump-2 reading
was a naive-unwrap artifact from the `s=−λ` crossing colliding with the on-line `s=0` tan-zero;
corrected.)* This is **not** a no-eigenvalue certificate, for three reasons:
1. **The symbol vanishes on the Plancherel line:** `R(0)=0` on `Re s=0` for *every* λ (`tan 0=0`), so
   the operator is **non-Fredholm** there — `Re s=0` is in the essential spectrum (the same fact as
   `edge_layer_estimate.md`: the essential spectrum sits on `Re=0`, and only an edge-unbounded weight
   could move it — which re-triggers cap_scaffold). The index is contour-side-ambiguous by 1.
2. **`W=3/2` is the raw symbol winding, not the index** (the `(s+λ)(s+2)` piece is the entire Γ-growth
   factor, not a bounded Toeplitz symbol; the tan-half is a jump-at-infinity). "Index 0" is uncomputed.
3. **Even index 0 ⊅ empty kernel.** `index = dim ker − dim coker`; the upgrade theorem is **Coburn's
   lemma** (`min(dim ker, dim ker*)=0` for a symbol nonzero on the contour) — **blocked** here by
   `R(0)=0`. Equivalent analytic form: prove the edge-recessive and far-field-recessive **Pincherle
   (minimal) solutions never coincide for `Re λ>0`**.

**Circularity (clean for the winding; firewall risk is elsewhere).** `W(λ)` is genuinely
spectrum-independent — it is the tautology `Re(−λ)<0 ⟺ Re λ>0` (location of the indicial root
`s=−λ`), computed from the explicit symbol. *Because* it is a tautology it carries no
no-eigenvalue content; the residual no-kernel step still routes through the `L²_W ⟺ B=0` firewall
until `τ(s)` is built. The winding also lives on the **symbol/recurrence variable**, dual-but-distinct
from the **eigenfunction-Mellin variable** (Taylor-coeff vs Mellin-transform, caveat C1) — a separate
bridge.

**Net (status UNCHANGED): `reduced-to-sharp-WH-factorization`.** The genuine spectral content lives
entirely in **`τ(s)`**. Sharpened remaining object:
> Construct the tan-Barnes `τ(s)` (`τ(s+1)/τ(s)=tan(πs/2)`) holomorphic and `L²`-decaying on the
> correct Mellin-inversion contour, with (1) the **side-choice** of the on-line zero `R(0)=0` resolved
> (it fixes the index at the `±1` level), and (2) a **Coburn-substitute / Pincherle-non-coincidence**
> step converting the index into `dim ker = 0` for `Re λ≥0` mod `𝒩`.

Gate: `validate/edge_wh_structure.py`.

## 5.3 The tan-Barnes `τ(s)` is ELEMENTARY — the construction (verified)

The remaining object (§5.2) was "construct the tan-Barnes `τ(s)`, `τ(s+1)/τ(s)=tan(πs/2)`,
holomorphic + `L²`-decaying on the inversion line." It is **elementary** (machine-verified to ~1e-31,
`validate/edge_tau_construction.py`):
`  τ(s) = π/sin(πs/2) = Γ(s/2)·Γ(1−s/2)`   (reflection),
since `τ(s+1) = π/sin(πs/2+π/2) = π/cos(πs/2)`, so `τ(s+1)/τ(s) = sin/cos = tan(πs/2)`. Hence the
eigenfunction's Mellin transform is an explicit **product of four ordinary Gammas**:
`  Ĝ_p(s) = Γ(s+λ)·Γ(s+2)·Γ(s/2)·Γ(1−s/2)`,   `Ĝ_p(s+1)/Ĝ_p(s) = (s+λ)(s+2)tan(πs/2) = R(s)`
(verified, ratio−R ~ 1e-30 at random complex `s,λ`), and `G(ξ)=(1/2πi)∫_C Ĝ_p(s)ξ^{−s}ds` is an
explicit **Fox‑H function**. The "research-grade" transcendental was classical all along.

**The Wiener–Hopf split is read straight off the four Gammas** (which Γ carries which poles):
`Γ(s/2)` → left poles `s=0,−2,−4,…`; `Γ(1−s/2)` → right poles `s=2,4,6,…`; `Γ(s+λ)` → `s=−λ−k`;
`Γ(s+2)` → `s=−2−k`. On the inversion line `Re s=0`, `Ĝ_p` decays `~e^{−3π|τ|/2}` (two Γ's at
`e^{−π|τ|/2}` + two half-Γ's at `e^{−π|τ|/4}`) — so it is `L²` on the line — with a single pole at
`s=0` (from `Γ(s/2)`) on the line = the marginal mode in `𝒩`.

**The point-spectrum question is now an explicit pole-killing-vs-`L²`-decay BALANCE** (general
solution `Ĝ=Ĝ_p·P`, `P` 1-periodic; eigenfunction needs `G` decaying at both ends): the edge `ξ→0`
is governed by the right poles (leading `ξ^{−2}` from `s=2` ⟹ edge-divergent unless killed), the far
field `ξ→∞` by the left poles (leading `ξ^{λ}` from `s=−λ` ⟹ divergent for `Re λ>0` unless killed).
Killing a residue class with a `sin`-type `P` costs `e^{π|τ|}` of growth against the `e^{−3π|τ|/2}`
decay; killing both bad classes costs `~e^{2π|τ|} > e^{3π|τ|/2}` ⟹ **no `L²` eigenfunction** — except
where poles **collide** (`λ∈ℤ`: `Γ(s+λ)` pole at `s=−λ` meets the `Γ(s/2)`/`Γ(s+2)` poles), which are
exactly the `𝒩`/symmetry/resonance modes. **Candidate conclusion: no `Re λ≥0` point eigenvalue mod
`𝒩` ⟹ the profile is stable.** ⚠️ The balance/growth-count is currently **heuristic** — under
adversarial determination (workflow `wf_57449c86`, residue + Pincherle + direct numerical `L²` scan +
circularity audit); status below pending that verdict. Gate: `validate/edge_tau_construction.py`.

**⚠️ CORRECTION (2026-06-14, adversarial verification): the BALANCE paragraph above has TWO bookkeeping
errors; the conclusion survives via an *exact integer* obstruction (§5.4).** (1) The edge↔far-field
pole assignment is **reversed**. Standard Mellin (`|ξ^{−s}|=ξ^{−Re s}`, `ξ=(Y−1)^{1/2}`) ⟹ far field
`ξ→∞` is governed by the **right** poles — leading `ξ^{−2}` from `s=2` (`Γ(1−s/2)`), which **decays**;
edge `ξ→0` by the **left** poles — leading the `s=0` const (`Γ(s/2)`) = the marginal `𝒩`-mode (raw
`Ĝ_p` decays at far field `ξ^{−2}` for *every* `λ`; numerically edge-exp `≈0`, far-exp `≈−2`,
`validate/edge_pointspectrum_rigid.py`). (2) The `sin(πs)`-type killer is **anti-periodic**
(`sin(π(s+1))=−sin(πs)`) so it does **not** preserve `Ĝ(s+1)=R(s)Ĝ(s)`; valid `P` are 1-periodic =
Laurent in `w=e^{2πis}` with **asymmetric** two-sided growth `e^{2πk_max|τ|}` (`τ→−∞`),
`e^{2π|k_min||τ|}` (`τ→+∞`). The corrected count is *sharper*, not weaker — see §5.4.

**Independent reproduction (2026-06-14, fresh mpmath, dps=50, no shared code).** All five construction
claims confirmed: (a) `τ(s+1)/τ(s)−tan(πs/2) ~ 9e-51`; (b) reflection `~1e-50`; (c)
`Ĝ_p(s+1)/Ĝ_p(s)−R(s) ~ 1e-49` over six `λ` (incl. complex) × twelve `s`; (d) the line-decay exponent
measured by a large-`τ` difference quotient (which cancels the polynomial Stirling prefactor) is
`−4.712330` vs `−3π/2 = −4.712389`, rel. err `1.3e-5`, converging monotonically — and `∫|Ĝ_p|²dτ`
finite once the `s=0` pole is excised, so genuinely `L²` on the line; the `s=0` residue is `2Γ(λ)`
(contour-verified to `~1e-48`). (e) The WH split is exact (`rgamma` pole-marking: `Γ(s/2)` left
`{0,−2,−4,…}`, `Γ(1−s/2)` right `{2,4,6,…}`, cleanly separated by `Re s=0`); the Fox-H integrand
reconstructed from `H^{1,3}_{3,1}[1/ξ | (1−λ,1),(−1,1),(1,½); (1,½)]` is **bitwise identical** to `Ĝ_p`
(diff `0.0`), and the edge `ξ^{−2}` term is the `s=2` residue verified by a signed two-contour
Mellin-Barnes inversion (`G_{c=3}−G_{c=1} = +Res_{s=2}[Ĝ_pξ^{−s}] = −18.536·ξ^{−2}`, agreement
`3.8e-34`). **Fox-H, not Meijer-G** (two distinct B-scales 1 and ½ — the `π/sin(πs/2)` half-scale block
is irreducible; only the `Γ(s+λ)Γ(s+2)` Pochhammer block is Meijer-reducible). No error found in (a)–(e).
The BALANCE growth-count is also arithmetically self-consistent (`e^{−3π|τ|/2}·e^{2π|τ|}=e^{+π|τ|/2}`
grows ⇒ not `L²`; one-sided `e^{−π|τ|/2}` still decays) — but remains heuristic for the reason the doc
states: a period-1 `sin(π(s−s₀))` kills an entire `ℤ`-spaced class at once, and that the *minimal*
multiplier killing exactly the two bad classes (edge spacing 2, far-field spacing 1) without re-creating
poles or destroying the `s=0` marginal mode costs *exactly* `e^{2π|τ|}` is the unproven
Pincherle/minimal-solution step. Gate (independent): `validate/edge_tau_independent_verify.py`.

## 5.4 The point-spectrum step, made rigorous — EXACT integer obstruction (2026-06-14)

The §5.3 BALANCE was heuristic *and* carried the two errors flagged above. The corrected, exact
version (gate `validate/edge_pointspectrum_rigid.py`, all claims direct mpmath, raw magnitudes /
`L²` integrals — no slope-sign shortcuts, which is where the heuristic hid its sign error):

**The exact homogeneous freedom.** `Ĝ=Ĝ_p·P` with `P(s+1)=P(s)` *exactly* — a 1-periodic
meromorphic factor, i.e. a function of `w=e^{2πis}`. (`sin(πs)` is anti-periodic and **invalid**: it
solves `Ĝ(s+1)=−R(s)Ĝ(s)`.) A 1-periodic `P=Σc_k w^k` (Laurent, extreme powers `k_min≤0≤k_max`) has,
on the line `s=c+iτ`:
- `τ→−∞` (`w→∞`): `|P|~e^{2πk_max|τ|}`,   `τ→+∞` (`w→0`): `|P|~e^{2π|k_min||τ|}`,
- and **exactly `d=k_max−k_min` zeros per period-strip** (roots of the degree-`d` Laurent polynomial)
  — so killing even one pole-coset needs `d≥1`, hence growth `≥e^{2π|τ|}` on at least one side.

**The two obligations.**
- *Line `L²` (Plancherel).* `Ĝ_p` decays **symmetrically** `~e^{−(3π/2)|τ|}` (λ-independent, reproduced
  to rel-err `~1e-5` across real/complex `λ`). `Ĝ_p·P∈L²(line)` ⟺ decays both sides ⟺
  `2πk_max<3π/2` **and** `2π|k_min|<3π/2` ⟺ `k_max<¾` **and** `|k_min|<¾` ⟺ `k_max=k_min=0` ⟺
  `d=0` ⟺ **`P` has no zeros.** (Verified: `∫|Ĝ_p·P|²` is `22.6` for `d=0` and astronomically
  divergent — `~10^{1.99×10^{34}}` — for every `d≥1` killer, symmetric or asymmetric.)
- *`L²_W` membership.* The `d=0` choice `P=1` leaves the edge `s=0` term = a nonzero **constant** at
  `ξ→0`. A const has **infinite** weighted edge mass: `∫_{−∞}^{z_cut}W·const²dz=+∞` for any
  bounded-at-edge `W` (`W→const`, the only admissible kind — `edge_layer_estimate.md`) over the
  half-infinite edge interval. So `L²_W` **requires** killing `s=0`, i.e. `d≥1`.

**The obstruction is `d≥1 ∧ d=0`, i.e. EMPTY.** No 1-periodic `P` makes `Ĝ_p·P` simultaneously line-`L²`
and `s=0`-killed. Hence **no `Re λ≥0` point eigenfunction** — *unless* an `Ĝ_p` pole **coincidence**
removes the `s=0` obligation. That happens exactly at `λ∈ℤ`: `Γ(s+λ)`'s tower merges onto the integer
cosets, the recurrence `a_{m+1}=(m+λ)(m+2)tan(πm/2)a_m` **terminates** (`tan(πm/2)=0` at even `m`),
`G` is a polynomial — these are precisely the `𝒩={marginal, dilation, symmetry}` modes that are
quotiented. So the only `Re λ≥0` candidates are `λ∈{0,1,2,…}=𝒩`; **modulo `𝒩`: empty ⟹ point-stable.**

**Why this beats the heuristic:** "`e^{2π}>e^{3π/2}`" is replaced by the exact integer inequality
`k_max<¾ ⟹ k_max=0` (no non-trivial integer in `[0,¾)`). The killer's growth is **rigid** — a
zero-free 1-periodic entire function is constant (`e^{as}` is not 1-periodic unless `a=0`), so there is
no theta-ratio / clever-`P` escape; the degree-vs-growth law is exact.

**Honest remaining gaps (this step is point-spectrum only, frozen generator):**
1. **Essential spectrum untouched.** `R(0)=0` on `Re s=0` for every `λ` (`tan0=0`) ⟹ non-Fredholm
   line ⟹ `Re=0` is in the essential spectrum. This point-mode argument says nothing about it; a
   **strict** spectral gap `Re σ≤−δ` still needs the edge-*unbounded* weight (re-triggers
   `cap_scaffold`), which is the separate open piece of `edge_layer_estimate.md`. So the earned claim
   is "**no `Re λ≥0` point eigenvalue mod `𝒩`**," not "spectral gap."
2. **The `L²_W`-edge-mass step** is the one conclusion-adjacent link; it is the honest content of the
   former firewall (`B=0 ⟺ L²_W`), now *computed* from the explicit `Ĝ_p` pole at `s=0` (the const
   has infinite half-line mass) rather than asserted — non-circular **iff** the `W→const` edge fact is
   used (it is). Worth a clean Lean/interval statement of "bounded-at-edge `W`, const ∉ `L²_W`."
3. **Frozen Ω-shape.** The nonlinear / co-moving-frame correction (`coevolve_*`) is a separate leg.

**Net status: `reduced-to-sharp-WH-factorization` → the point-spectrum half is now an EXACT integer
obstruction** (modulo gap 2's clean weighted-norm statement); the essential-spectrum gap (1) is the
genuine remaining analytic object for a *strict* stability margin. Gate:
`validate/edge_pointspectrum_rigid.py`.

## 5.5 Route 1 (residue / pole-collision): the BALANCE heuristic is sharpened to a structural bound

A fresh adversarial pass (2026-06-14, independent code `validate/edge_pole_collision.py`, dps=30)
**confirmed every prior claim** and **closed two loopholes** in the §5.3 BALANCE heuristic, turning it
from a growth-count slogan into a near-rigorous structural argument. The remaining gap is now a single
named lemma, not a hand-wave.

**Edge residue (exact, analytic).** `Res_{s=2} Ĝ_p = −2·Γ(2+λ)Γ(4)Γ(1) = −12 Γ(2+λ)`, **nonzero for
every generic `λ`** (numeric = analytic to `~1e-30`; `λ=0.7 ⟹ −18.536`, matching the §5.3 reading). So
with `P=1` the edge term `ξ^{−2}` genuinely diverges (`ξ→0`). With `P=1` the eigenfunction does not exist.

**Pole orders (the heuristic was OPTIMISTIC).** Right poles `s=2,4,…` SIMPLE (class `{0}` mod 1); left
`s=0` simple (marginal, in `𝒩`); **`s=−2,−4,…` are DOUBLE** (`Γ(s/2)`&`Γ(s+2)`, class `{0}`) — TRIPLE at
integer `λ`; `s=−λ−k` simple (class `{−λ}` mod 1). The far-field leading divergence is therefore not
`ξ^λ` (from `s=−λ`) but the **stronger `ξ^2 log ξ`** from the `s=−2` double pole — this makes the
far-field worse, so the conclusion is more robust, not less.

**Loophole 1 closed — `P` must be ENTIRE (structural theorem).** A 1-periodic meromorphic `P` has ALL
its poles/zeros with `Re s ∈ (−½,½]` (one period). `Ĝ_p` is **zero-free** (a product of Γ's, which never
vanish). Hence any pole of `P` lands on/beside the inversion line `Re s=0` and **survives** in `Ĝ_p·P`
(`Re s_p = arg(w_p)/2π`, always in `(−½,½]` — verified: a real-positive `w_p` pole sits exactly at
`Re s=0`). So no meromorphic `P` can use off-strip poles to cut its line-growth: **`P` must be an entire
1-periodic function** `P(s)=Σ_n c_n e^{2πi n s}` (a Laurent polynomial in `w=e^{2πi s}`). This is the step
the §5.3 heuristic silently assumed.

**Loophole 2 closed — per-zero growth is exactly `π|τ|` (centered Fourier).** An entire 1-periodic `P`
with `D` zeros per period (with multiplicity) has minimal SYMMETRIC two-sided line-growth `e^{πD|τ|}`
(Fourier support `n∈[−D/2,D/2]`; verified: `sin(πs)` has `D=1`, slope `+π=+3.14159`).

**The inequality (arithmetic on two spectrum-free exponents).** `Ĝ_p` line-decay `= −3π/2` (Stirling).
Required killer degree: integer class `{0}` needs order `max(edge 1, far-double 2)=2`; class `{−λ}` needs
order 1 ⟹ `D_req ≥ 3` (generic `λ`). Line-slope of `Ĝ_p·P = −3π/2 + πD_req ≥ −3π/2 + 2π = +π/2 > 0`
(even the over-optimistic `D=2` ignoring the far double pole loses; `D=3` gives `+3π/2`). **Positive slope
⟹ not `L²` ⟹ no eigenfunction, for every non-integer `Re λ≥0`** (numerically confirmed: `D=2` slope
`≈+1.6`, `D=3` slope `≈+4.77`, across `Re λ∈[0.05,3.6]` incl. complex `0.5+1.5i`). A single-class (`D=1`)
killer is `L²`-on-line but **insufficient**: `Ĝ_p·sin(πs)` leaves the `s=−λ` far-field pole alive
(`G(ξ)~ξ^λ` at `ξ→∞`) — the eigenfunction genuinely needs BOTH ends, forcing `D≥2`.

**Collision = exactly `𝒩`.** Integer `λ` merges `{−λ}` into `{0}` (one class) but the left orders go UP
(triple) — these modes are `λ=1` (time-translation symmetry, removed by `𝒩`) and `λ` even (`tan(πζ)=0` ⟹
REGULAR edge modes `λ≤0`). No `Re λ>0` point eigenvalue survives mod `𝒩`.

**Firewall (clean).** The two exponents `−3π/2` (Stirling) and `+π`/zero (periodicity) are computed
independently of any stability assumption; the conclusion is the **arithmetic** `D_req·π ≥ 2π > 3π/2`. No
`L²_W ⟹ B=0` substitution is used; the "`P` entire" step comes from PERIODICITY, not from benign-spectrum.
Audit gate `validate/edge_synthesis_firewall_audit.py`: the only conclusion-adjacent link
(`const ∉ L²_W`) is re-expressed as "`L²_W` (= Mellin-`L²` on `Re s=0`) ⟹ no pole at `s=0`" — a
*property of the weight* computed from `Ĝ_p`'s explicit `s=0` residue `2Γ(λ)`, NOT injected from the
benign spectrum. The content is the **rigidity**: the 4-Gamma product cannot be simultaneously
(i) recurrence-correct, (ii) Mellin-`L²` on the line, (iii) `s=0`-pole-free.

**Remaining gap (sharply stated).** The last unproven step is the **Pincherle/minimal-solution lemma**:
that no entire 1-periodic `P` — possibly with EXTRA zeros arranged to lower one-sided growth — can beat
`e^{πD|τ|}` for genuine *two-sided* `L²` control. The centering argument bounds this and the `D≥2` margin
makes the verdict robust, but until that lemma is a stated-and-proved proposition the status stays
`reduced-to-sharp-WH-factorization` (now: **reduced-to-one-Pincherle-lemma**). Gate:
`validate/edge_pole_collision.py`.

**Candidate conclusion (earned confidence: strong-numeric + structural, NOT proof).** No `Re λ≥0` point
eigenvalue mod `𝒩` ⟹ the Hou–Luo profile is spectrally stable at the edge — pending the one Pincherle lemma.

## 5.5 The Pincherle lemma is CLOSED by Liouville (2026-06-14, Route 3 independent verification)

A fresh adversarial pass (independent code, no shared helpers — gates `validate/edge_route3_L2_scan.py`,
`edge_route3_minimal_multiplier.py`, `edge_route3_critical_line.py`, `edge_route3_periodic_correct.py`;
dps 18–30) **independently re-derived the entire corrected chain** (both the reversed edge/far assignment
*and* the `sin(πs)`-anti-periodicity error were re-found from scratch, then matched §5.4) and **closes the
one named gap** — the "Pincherle/minimal-solution lemma … possibly with EXTRA zeros arranged to lower
one-sided growth" left open in the Route-1 §5.4.

**The lemma, now proved (Liouville).** Let `P` be the homogeneous freedom: `P(s+1)=P(s)` and `P` entire
(loophole 1, §5.4: a meromorphic `P` keeps its poles in the period-strip on the line, since `Ĝ_p` is
zero-free — so `P` MUST be entire). Then `P(s)=F(w)`, `w=e^{2πis}`, with `F` **entire** in `w∈ℂ∖{0}`.
On the inversion line `s=σ+iτ`, `w=e^{2πiσ}e^{−2πτ}` traces the punctured ray `|w|: 0←(τ→+∞) … (τ→−∞)→∞`.
So:
> `P` bounded on **both** sides of the line ⟺ `F` bounded as `|w|→0` **and** `|w|→∞` ⟺ `F` extends to a
> bounded entire function on all of `ℂ` ⟺ (Liouville) `F ≡ const` ⟺ `P` zero-free constant (`d=0`).

Hence **every nonconstant entire 1-periodic `P` grows at least like `|w|^{1}=e^{2π|τ|}` on at least one
side** — NO arrangement of extra zeros escapes (that was the open loophole). The minimal one-sided growth
per zero is therefore the integer rate `2π`, not the half-integer-Fourier `π` the centering heuristic
suggested (verified: a single-zero integer-Laurent `e^{2πis}−1` has max one-sided slope `+2π=6.283`, while
the symmetric `w+w^{−1}−c` (`d=2`) grows `e^{2π|τ|}` on BOTH sides — `validate/edge_route3_periodic_correct.py`).

**Consequence — the margin holds already at `d=1`.** `Ĝ_p` decays symmetrically `e^{−3π/2|τ|}`; one zero
costs `e^{+2π|τ|}` one-sided ⟹ net `e^{+π/2|τ|}` GROWTH ⟹ not line-`L²`. So even a *single*-coset killer
fails on the line, before the far-field double-pole counting — the `2π > 3π/2` integer inequality
(`k_max<¾ ⟹ k_max=0`) is the whole story, and it is now backed by Liouville rather than a centered-Fourier
bound. The two obligations remain exclusive: `s=0`-kill needs `d≥1` (else infinite weighted edge mass);
line-`L²` needs `d=0`. **Empty ⟹ no `Re λ≥0` point eigenfunction mod `𝒩`.** Verified end-to-end: `sin(πs)`
gives ratio `−R` (inadmissible); `Ĝ_p·sin` is line-`L²` but leaves `ξ^λ` far-blow-up (a true eigenfunction
needs both ends); `Ĝ_p·sin·sin` and every entire `d≥1` killer grow on the line.

**Status upgrade: `reduced-to-one-Pincherle-lemma` → point-spectrum half CLOSED (structural proof).** The
Liouville step is elementary and spectrum-independent (firewall clean — it uses only periodicity + the
`Ĝ_p` zero-freeness + Stirling decay, never a stability assumption). The THREE honest gaps from §5.4
**remain and are unchanged**: (1) the **essential** spectrum on `Re=0` (`R(0)=0`, non-Fredholm line) is
untouched — a *strict* gap `Re σ≤−δ` still needs the edge-unbounded weight (`cap_scaffold`); (2) the
`const ∉ L²_W` edge-mass step deserves a clean interval/Lean statement (`W→const` at the edge over a
half-infinite `z`); (3) frozen Ω-shape, nonlinear correction separate. So the earned claim is precisely
**"no `Re λ≥0` POINT eigenvalue modulo `𝒩`"** — point-spectral stability of the frozen generator — not a
strict spectral gap. THR-01 stays `speculative` on the *strict* margin; the point-spectrum obstruction is
now a structural argument with no remaining hand-wave. Gate: `validate/edge_route3_critical_line.py`,
`validate/edge_route3_periodic_correct.py`.

## 5.6 Synthesis-judge verification (2026-06-14, `wf` synthesis pass) — point-spectrum half CONFIRMED, one inter-route discrepancy resolved

A fresh independent reproduction (all gates re-run + five new adversarial probes
`validate/edge_synthesis_*.py`, dps 25–30, no shared helpers) **confirms the construction exactly
and the corrected point-spectrum obstruction structurally**, and resolves a numerical discrepancy
*between* the two prior routes that did not affect the verdict.

**Construction (a)–(e): reproduced to `~1e-49…1e-53`** — `τ` elementary, `Ĝ_p` 4-Gamma, recurrence,
symmetric line-decay `3π/2` (rel-err `~1e-4` to `1e-5` across real/complex/negative `λ`), `s=0`
residue `2Γ(λ)` (to `~1e-38`), WH split, Fox-H bitwise-identical, edge `ξ^{−2}=Res_{s=2}` to `1e-53`.
No error in the breakthrough construction.

**The π-vs-2π per-zero discrepancy is RESOLVED in favor of 2π (Route 3 §5.5 correct; Route 1 §5.4
overstated the cheapness).** Route 1 (`edge_pole_collision.py`, step 4) claimed per-zero *symmetric*
growth `e^{π|τ|}` (centered Fourier, citing `sin(πs)` `D=1 → π`); Route 3 (`edge_route3_periodic_correct.py`)
and the rigid gate claimed *one-sided* `e^{2π|τ|}` (integer-Laurent). Direct measurement
(`edge_synthesis_per_zero_growth.py`) settles it: the only object giving symmetric `π` at one zero is
`sin(πs)`, which is **anti-periodic and inadmissible** (it solves `Ĝ(s+1)=−R(s)Ĝ(s)`, ratio/R `=−1`
exactly). Every *admissible* (genuinely 1-periodic, `=` Laurent in `w=e^{2πis}`) single-zero killer —
`w−1`, `w−c`, `1−w^{−1}` — has growing-side rate **exactly `2π`** and decays the other side. So
Route 1's inequality `−3π/2 + πD_req` used an inadmissible per-zero constant; the correct
budget is the integer one-sided rate `2π`. **The verdict is unaffected** — Route 1 still lands at
`D=3` (positive slope), Route 3 reaches the same conclusion already at `d=1` — but Route 3's
`2π > 3π/2` form is the correct and sharper one. (Route 1's §5.4 §"GROWTH PER ZERO = exactly `π|τ|`"
line should be read as superseded by §5.5's `2π`.)

**Pincherle "extra-zeros" loophole CLOSED, independently (`edge_synthesis_pincherle_loophole.py`).**
An entire 1-periodic `P` is a *Laurent polynomial* in `w` (else essential singularity at `w=0` or
`w=∞`), so any `P` with ≥1 zero has growing-side Laurent degree ≥1 ⟹ rate ≥`2π` on at least one side;
no arrangement of extra zeros lowers the binding-side rate below `2π` (verified: `w−1`, `w−0.5`, `w−5`,
`1−w^{−1}` all measure `2π`; the would-be `π` object `(w−1)/w^{1/2}` has a nonzero 1-periodicity
residual `~1e3` ⟹ not 1-periodic). This is the Liouville closure of §5.5, re-derived from scratch.

**Firewall audit (`edge_synthesis_firewall_audit.py` + `edge_synthesis_admissible_killer.py`): clean,
with the honest caveat made explicit.** The logical skeleton is: recurrence-correct ⟹
`Ĝ=Ĝ_p·P`, `P` 1-periodic entire; `L²_W` membership `≡` (ii) line-`L²` on `Re s=0` **and**
(iii) holomorphic on the line (no `s=0` pole). The 4-Gamma rigidity forces (ii)`⟹d=0` but
(iii)`⟹d≥1` — **empty intersection** (exhaustive integer-Laurent scan to deg 4 finds no `P`
satisfying both; the admissible `s=0`-kill `w−1` alone breaks line-`L²` by `2π−3π/2=π/2` one-sided,
so the far-field `ξ^λ` step in `edge_route3_critical_line.py` is *redundant* — and that script's
`sin(πs)` "candidate eigenfunction" is logically void since `sin` is inadmissible, though it reaches
the same place). The conclusion-adjacent fact "forced edge constant `B=2Γ(λ)≠0 ∉ L²_W`" is **computed**
from the explicit `s=0` residue and the Mellin-Plancherel weight, *not* injected from benign spectrum
— non-circular. **Honest framing:** items (ii)+(iii) *are* the `L²_W`-membership definition, so the
genuine content is the **rigidity** that the explicit 4-Gamma product cannot meet both simultaneously;
that rigidity is firewall-clean (Liouville + Laurent + `Ĝ_p` zero-freeness + Stirling, no stability
assumption). The threshold `3π/2 < 2π` (`k<3/4 ⟹ k=0`, no positive integer below `3/4`) is **exact
with zero slack** (`edge_synthesis_threshold.py`).

**Synthesis verdict.** Construction: **sound, fully reproduced, machine-exact.** Point-spectrum
conclusion ("no `Re λ≥0` eigenvalue mod `𝒩`"): the three routes (residue/collision, Pincherle/Liouville,
direct numerical `L²` scan) **agree**, and after the π→2π correction the structural argument has **no
remaining hand-wave for the point spectrum** — the rigidity is an exact integer obstruction backed by
Liouville. It is **not** a hand-wavy growth count and does **not** re-encode the firewall (the firewall
step is the *definitional* `const ∉ L²_W`, computed from the residue). The honest residual is what §5.4–5.5
already state: **(1)** the **essential** spectrum on `Re=0` (`R(0)=0`, non-Fredholm) is untouched — a
*strict* gap `Re σ≤−δ` is still open and needs the edge-unbounded weight; **(2)** the `const ∉ L²_W`
edge fact deserves a clean Lean/interval statement; **(3)** frozen Ω-shape, nonlinear leg separate.
**Earned status: point-spectrum half is a structural argument (strong-numeric + structural, three routes
agree); the strict spectral gap remains open. THR-01 stays `speculative` on the strict margin.** Gates:
`validate/edge_synthesis_{per_zero_growth,pincherle_loophole,firewall_audit,admissible_killer,threshold}.py`.

## 6. Circularity firewall (holds; read before any write-up)

The symbol derivation (§1–§4) is **non-circular**: an explicit meromorphic symbol computed from
the operator definition, no stability injected (8/8 verdicts `circular=false`). The single
conclusion-adjacent step is `B = H(g)(1)/(2λ) = 0` on `Re λ≥0` modes — and it is **identical to
`L²_W` membership** (a nonzero forced constant `∉ L²_W` on the infinite `z`-line), NOT an
independent scalar lemma. The route stays non-circular **iff `B = a(0)` is genuinely computed
from the Wiener–Hopf product** (an independent object). It re-enters the NS-C6/SQG trap the moment
a write-up substitutes "the discrete spectrum is benign / the mode is admissible ⟹ `B=0`" for the
actual index computation. The continuum anchor `H[(Y−1)^{−1/2}]=0` is exactly what the **discrete**
Hilbert violates at the edge node (`−0.445·e^{−z_min/2}`), so this step provably **cannot** be
done numerically.

## 7. Numerics reconciliation (all consistent — OBS-03, OBS-04)

`λ=−2ζ` ladder; `+40 = ζ=½` resonance (κ-crushable, g-unbounded branch); the marginal transport
essential spectrum on `Re=0` (`max Re → 0` as `N→∞`, `Im ~` Nyquist, weight-invariant on a finite
grid) — a direct eigenvalue computation cannot exhibit the gap (it needs the edge-divergent
weight). The `+1.7` "instability" (edge_bc off) is a `z_min`-drifting boundary artifact (the
discretized inadmissible homogeneous branch). None contradicts the symbol picture.

## 8. Status and next step

**Status: `reduced-to-sharp-WH-factorization`.** From "prove no `Re λ≥0` spectrum of a strongly
non-normal nonlocal operator (energy methods provably dead)" down to **one** bounded analytic
object — the Wiener–Hopf index of a meromorphic, single-valued Mellin difference equation. THR-01
stays `speculative`; the index is named, not evaluated.

**Next concrete step.** (i) Argument-principle index in-strip: fix `Re λ ∈ [0,1)`, truncate the
down-ladder product (the `1/n!`-type decay converges fast), evaluate `(1/2πi)∮ d log a(s)` around
the strip — verify `0` (no eigenvalue) where the symbol is regular. (ii) Push `Re λ ≥ 1`: recompute
the index across the strip boundary via the analytic continuation (Plemelj/Gakhov for `cot`-type
symbols on a graded lattice) — the one step that needs the continuation, not the integral.
**Honest scoping:** scalar WH factorization of a `cot`/`tan` symbol on a single line is classical
(Gohberg–Krein); on the **half-integer-graded lattice with two matched strips and the `𝒩`-quotient**
it is most likely a **small research problem**, not an off-the-shelf citation. That graded-lattice
factorization is the leap.
