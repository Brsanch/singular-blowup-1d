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
