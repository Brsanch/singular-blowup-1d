# desk/ attack doc — the edge-layer estimate (LEM-06 / the (★) leap)

Date 2026-06-13. OPEN-GAP desk doc (proof-obligation-ledger convention) for the
edge-layer leg of (★)/LEM-06 — the last analytic piece of the Conjecture-2.4
stability proof for the singular Hou–Luo profile. Method: 6 independent derivation
lenses, each adversarially verified by 3 skeptics (sign/constant, circularity,
numerics), plus an independent numerical reproduction (`validate/edge_layer_obstruction.py`).

**Bottom line (honest, and it is a NEGATIVE result): the edge layer does NOT close by
any weighted-energy / dissipativity / Schur / flux method, and this is now *proven*,
not merely "not yet done". The leap is no longer an energy estimate — it is a
non-self-adjoint spectral/Fredholm bound. Status: `open-obstruction`.**

---

## 0. The reduction (recap, established)

Frozen Ω-shape generator on `z = ln(Y-1) ∈ (−∞,∞)`, `Y = 1+e^z`:
`M_GG = −2∂_z + F`, `F = ½Q_op`, with the CAUSAL exponential-averaging forcing
`(Q_op G)(z) = ∫_{−∞}^z e^{−(z−z')} H[e^{−·/2}G](z') dz'`, `H` = line Hilbert on `Y>1`.
Weighted energy `‖G‖²_W = ∫ W G² dz`, anisotropic `W = e^{−φ}`, `φ'(z) = 2γ·½(1+tanh((z−z0)/s))`
→ 0 at the edge, → 2γ in the bulk. Transport is exactly dissipative (Lean: `weighted_ibp`):
`⟨G,−2∂_z G⟩_W = −∫φ'WG² − [WG²]`, but `φ'→0` at the edge so the local edge dissipation → 0.

Split at `z_cut`: the **bulk** (`z>z_cut`) coercivity is numerically solid and
**depth-stable** (`A_bulk = −0.41`(z>−6)/`−0.65`(z>−4), constant across `z_min=−10…−16`;
reproduced). The **edge layer** (`z<z_cut`) is the open piece, and is what this doc settles.

## 1. Main result of this attack — the energy route is dead at the edge

All six lenses converge (18/18 adversarial verdicts non-circular): the edge-layer
quadratic form `⟨G,FG⟩_{W,z<z_cut}` is **not** dominated by the available dissipation in
**any** bounded-at-edge weight. The mechanism is a single hard scaling fact:

> the nonlocal forcing density inherits the base-profile amplitude `(Y−1)^{−1/2} = e^{−z/2}`
> at the edge. A bounded-at-edge weight (`W→const`, the only kind that doesn't blow up the
> edge↔bulk Hilbert coupling) cannot cancel an `e^{−z/2}` density. The only weight that
> would dominate it has `φ'≳1` (≳`e^{−z}` in `W`) at the edge — i.e. the **convective weight
> already ruled out** in `findings_cap_scaffold.md` (its field-of-values diverges `C_force=10→29→83`).
> The two failure modes are **dual**: pushing dissipation edgeward to beat the forcing
> re-triggers the convective non-commutation. No admissible window exists.

Quantitatively (lens-derived, then **independently reproduced** here):

- **Self-similar collapse.** The full symmetrized energy abscissa diverges geometrically
  with edge depth at *exactly* the profile rate: `A_full · e^{z_min/2} = 0.099 ± 0.001`
  across `z_min = −10/−12/−14/−16` (`A_full = +14.7/+40.2/+108/+291`). This is a **continuum**
  exponent `(Y−1)^{−1/2}`, not a grid rate (a grid rate would scale with `N`, not `z_min`).
- **Schur / flux / FoV all share this divergence.** kernel-schur: edge Schur norm `∝ e^{−z_cut/2}`
  (`178→504→1381→3757→10214`). flux-matching: the `z_cut` boundary flux is **trace-bounded**
  (`C_trace = 0.96`, O(1), converged) while the forcing range is O(40) and diverging — ratio→0,
  the flux mechanism is real, correctly signed, and **quantitatively negligible** (`|flux|/|forcing|~10⁻⁵`).
  causal-volterra: the finite-rank `½c₁` functional is **unbounded** on the slow-decay class.

## 2. Numerical reconciliation — artifact vs. real (the dispute the skeptics surfaced)

The lenses split on whether `+40` is a grid artifact or a continuum obstruction. It is
**both, and neither is the whole story** — resolved by `validate/edge_layer_obstruction.py`:

| test | result | reading |
|---|---|---|
| `A_full·e^{z_min/2}` | `0.099` const, `z_min=−10…−16` | **continuum** `(Y−1)^{−1/2}` resonance, not a grid rate |
| bulk `A(z>z_cut)` | `−0.41`/`−0.65`, depth-independent | divergence is **localized to the edge layer**; bulk is genuinely coercive |
| g-mass penalty `−κ·diag(e^{−z})` | `κ=10⁻³ → A_full = −0.24` (101% crush) | the `+40` lives **entirely on the g-unbounded `(Y−1)^{−1/2}` branch** |
| `H(profile)[node0]` | `−0.445·e^{−z_min/2}` (continuum = 0) | the **discrete Hilbert mis-resolves** the profile at the deepest node — genuine contamination, same `e^{−z/2}` rate |
| N-drift at fixed `z_min` | `A·e^{−6}: 0.093→0.107` over N=900→2800 | slow (~15%/3×), sub-log: the discrete echo riding on the continuum 0.099 |

**Synthesis:** the `+40` is a **continuum non-normal profile-resonance carried by the
inadmissible g-unbounded branch** (the indicial branch `Re λ ∈ [−1,0)` that is `L²_W`-admissible
but `g = e^{−z/2}G ∉ L²(dz)`). It is excluded by admissibility — but the repo enforces
admissibility only numerically (one dropped node / the κ-penalty here), **not by proof**.
On the genuinely admissible (g-bounded) class the abscissa is **negative but marginal**
(`≈ −10⁻³`, = the true spectral max `−8×10⁻⁴`), **not** a robust coercive `−δ`. So even with
admissibility imposed, the *energy* estimate yields only marginality, because the operator is
**strongly non-normal** at the edge (`A_FoV − spec ≈ +40`).

## 2.1 Spectral-side cross-check (2026-06-14, `validate/edge_spectrum_continuation.py`, OBS-04)

The cheap discriminating run (track the rightmost `M_GG` spectrum as `z_min→−18`, `N→2800`)
both supports stability and sharpens *why the numerics cannot close it*:

- **No `Re λ > 0` eigenvalue** appears anywhere in the window (`max Re ≈ −10⁻³`, `0` modes at
  `Re>0`), `→ 0` as `N→∞`. Positive evidence for stability.
- **But that maximum is the MARGINAL transport essential spectrum on `Re=0`, not a gap.** It is
  κ-stable (an edge g-mass damping `−κ·diag(e^{−z})` that crushes the `+40` abscissa to negative
  at `κ=10⁻³`, OBS-03, does **not** move it even at `κ=10⁻²`), delocalized, conjugate-paired,
  `→0` with `N`. So the repo-wide "benign spectrum `≈−10⁻³`" is the **marginal essential
  spectrum**, *not* a spectral gap — and the `+40` field-of-values (κ-crushable, inadmissible
  branch) and this essential spectrum are **distinct objects**.
- **Consequence:** the asymptotic-stability gap (`Re spec ≤ −δ`) is irreducibly a *weighted-space*
  property — only a weight shifts the `Re=0` essential spectrum left, and that weight is exactly
  what diverges at the edge. **A direct eigenvalue computation cannot exhibit the gap**, which
  confirms the leap is not a numerical question.
- **Honest null:** the membership question `B = H(g)(1)/(2λ) = 0` is **not** settled numerically.
  A smoothness (hi-k) filter does not robustly isolate the discrete point spectrum on this strongly
  non-normal operator (the "rightmost point mode" flips `−0.7 ↔ −22` with the scoring convention),
  and the edge-decay / `c₁` readouts are corrupted by oscillation + the discrete-Hilbert bias. The
  eigenvalue method inherits the same edge ambiguity as the energy method — reinforcing §5.

## 3. Route table (6 lenses)

| lens | requirement | what would close it | difficulty | verdict |
|---|---|---|---|---|
| kernel-schur | bounded Schur norm of `|K|` as `z_cut→−∞` | signed cancellation (signed/abs row sum ~4×) + `z_cut` flux identity — **outside Schur** | high | sound NEGATIVE; Schur discards the cancellation |
| causal-volterra | `c₁=H(g)(1)` bounded on admissible `H¹_W` | a-priori edge decay `a>½` (`Re λ<−1`); or `φ'≳1` at edge (tension w/ cap_scaffold) | high | sound; proves `c₁` **unbounded** on `a∈(0,½]` |
| matched-asymptotic | forced amplitude `B=0` by matching | Fredholm proof `H(g)(1)=0` on `Re λ≥0` modes | high | **partial, over-credited** — `B=0 ⟺ H(g)(1)=0` is the *same* condition as `L²_W` membership, not a new lemma (skeptics) |
| L2-membership | on `ker ℓ`, residual forcing < edge dissipation | — | high | **partial; "closes on ker ℓ" is FALSE** (dimensional error, wrong `φ'` exponent) — `ℓ=0` necessary, not sufficient |
| flux-matching | `z_cut` flux absorbs edge forcing | — | low | sound NEGATIVE: flux trace-bounded `0.96`, forcing diverges; flux is negligible |
| numerics-faithful | resolvent/pseudospectral bound on the non-normal edge block | verified `‖(L−λ)^{−1}‖≤M` (cap_scaffold route 2) **or** Mellin proof `B=0` | **leap** | best-calibrated; pins the leap, stops there |

## 4. The circularity firewall (this is the NS/SQG failure mode — read before any write-up)

`H(g)(1) = 0` on `Re λ ≥ 0` modes is **conclusion-adjacent**: it asserts no unstable mode
carries a forced edge singularity, i.e. "the edge is benign modulo the bulk" ≈ the stability
conclusion. Two facts keep the program honest:

1. **It is automatic, not an extra constraint.** A `Re λ≥0` eigenmode has its homogeneous edge
   branch `e^{−λz/2}` blowing up at the edge (inadmissible, dropped), leaving the forced constant
   `G→B`; but a nonzero constant `∉ L²_W` (infinite z-interval, `W→const`). So **`L²_W` membership
   alone forces `B=0`** ⟺ `H(g)(1)=0`. The matched-asymptotic "extra codim-1 constraint" framing
   is **over-credited** (skeptics; the paraphrase antipattern this codebase polices) — it re-expresses
   the `L²_W` domain condition, it does not reduce the problem to a new scalar lemma.
2. **Therefore the real question is the original global one**, now correctly scoped: *does the `L²_W`
   realization of `M_GG ⊖ 𝒩` have any `Re λ ≥ 0` eigenvalue?* — and every weighted-energy instrument
   to answer it provably **fails at the edge** (§1–§2). Asserting "no, because the spectrum is
   benign / the mode is admissible" **without an independent non-energy bound is exactly the NS-C6 /
   SQG-vacuous-hypothesis circular step.** No lens took it; a write-up must not either.

## 5. The sharpest obstruction (maximum resolution — the deliverable for the leap)

> **Prove the `L²_W` realization of `M_GG = −2∂_z + ½Q_op`, restricted to `H¹_W ⊖ 𝒩`
> (`𝒩 = span{ξ@λ=1, μ, ν}`), has no spectrum in `Re λ ≥ 0` — by a NON-ENERGY method**, since
> every bounded-at-edge weighted-energy estimate provably diverges at the edge (`A_FoV = 0.099·e^{−z_min/2}`,
> reproduced), and the convective weight that would dominate re-triggers the `cap_scaffold`
> non-commutation. Two concrete handles:
>
> **(a) Resolvent / pseudospectral bound** (cap_scaffold route 2): verified `‖(M_GG − λ)^{−1}‖_W ≤ M`
> on `Re λ ≥ 0` (mod `𝒩`). Non-normality-robust; must be in the weighted space (essential spectrum
> sits on `Re=0` unweighted) and so inherits the `e^{γz}` conditioning — interval arithmetic, not floats.
>
> **(b) Mellin / Wiener–Hopf symbol of `½Q_op` at the edge indicial exponent `r = −(λ+1)/2`.**
> Write the edge problem `(M_GG − λ)G = 0` as `z→−∞` as the inhomogeneous `G = A e^{−λz/2} + B`
> with the nonlocal source; the causal Volterra kernel `Q_op` composed with the line-Hilbert has a
> computable Mellin symbol whose value at `r` decides whether a nonzero forced constant `B` is
> admissible for `Re λ ≥ 0`. The load-bearing **continuum** fact is `H((Y−1)^{−1/2}) = 0` on the
> half-line (G–R 3.222.2) — the exact identity the **discrete** operator violates at the edge
> (`H(profile)[node0] = −0.445·e^{−z_min/2} ≠ 0`), so **the numerics actively mislead here and the
> step cannot be done numerically.**

This is non-interpolative because: (i) it is a solvability/spectral statement on a strongly
non-normal nonlocal (causal-Volterra ∘ line-Hilbert) operator — the entire energy toolkit is
*provably* inapplicable, not merely unattempted; (ii) it must use a continuum Hilbert identity the
discretization breaks; (iii) it cannot be justified from the benign spectrum without circularity.

## 6. Concrete next computations

1. **Discriminating numerics (cheap, settles the conjecture before the proof):** continue the
   rightmost `L²_W` eigenmodes as `z_min` deepens (`−12→−18`) **and** `N` refines, tracking the
   edge-decay rate `a` and `c₁ = H(g)(1)`. If `a` stays `≫ ½` and `c₁ → 0` under *both* limits,
   strong evidence the global coupling self-consistently enforces `B=0`. (Still numerical — the
   discrete Hilbert bias is the caveat; report it.)
2. **The Mellin symbol of `½Q_op` at `r = −(λ+1)/2`** — the one new analytic object. Its zeros in
   `Re λ ≥ 0` are the candidate instabilities; proving it has none (mod `𝒩`) closes the edge layer.
3. Reference architecture for the resolvent/CAP route: `taklab-org/RigorousIntegrator_higher_dimensional_PDEs`
   (Gershgorin-after-similarity, eigenpair enclosures) — but the transport-dominated, **non-sectorial**
   `M_GG` is outside their sectorial-semigroup tail bound; that gap is the new work (NEXT_SESSION §B-2b).

## 7. What Bryan must leap

The whole Conjecture-2.4 edge obligation reduces to a **single non-self-adjoint spectral statement**:
no `Re λ ≥ 0` spectrum of the `L²_W` realization `M_GG ⊖ 𝒩`, proven without a weighted-energy
inequality (all of which die at the edge) and without asserting `H(g)(1)=0` from the spectrum
(circular). The concrete handle is the Mellin/Wiener–Hopf symbol at `r = −(λ+1)/2`, using the
continuum identity `H((Y−1)^{−1/2}) = 0`. Everything else — the bulk coercivity (`δ=0.02…0.78`),
the transport IBP signs (Lean), the flux outflow sign, the field-of-values dead-end, the
profile-resonance rate `0.099·e^{−z_min/2}` — is established and *reduces the edge obligation to
this one object*.
