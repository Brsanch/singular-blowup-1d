# Findings — W3 by exact conjugation to damped transport

**Date:** 2026-11-21
**Probe:** `validate/edge_w3_conjugation_check.py` (WC-1..WC-3, exit 0, about 30 s, NumPy/SciPy).
**Review form:** independent — one Claude reviewer with fresh context, given the claim and told to break
it; own derivations and own numerics (scripts kept unedited in `validate/review_2026_11_21_w3/`; the
implementer has not re-run them). The implementer is a different Claude session. Not expert review.
Nothing here is machine-checked in Lean.
**Tier:** proved-on-paper (short argument, sections 1–2) + one independent adversarial pass that did
not break it. No ledger status is changed by this note.
**Does not claim:** anything about the Gauge-(a) operator `−2∂_z + ½Q̃_op` that the manuscript's
Theorem 2 evolves · that solutions satisfy the original system across `X = 1` (review item N1 stands)
· Conjecture 2.4 · blowup of the PDE.

---

## 0. Result in four lines

1. With `N(s) = cos(πs)Γ(2s+2)/π`, the symbol identity `N(s+½)c(s+½) = −N(s)` holds.
2. Hence on `L²_σ`, `σ ∈ (0,½)`: `M_σ = 𝒩 (−2∂_z − e^{−z/2}) 𝒩⁻¹`, where `𝒩` is the multiplier with
   symbol `N`, an isomorphism `H^p_σ → L²_σ`, `p = 2σ + 3/2`.
3. The conjugated equation is transport with edge damping and is solved in closed form, so `M_σ`
   generates a C₀-semigroup with `‖S(τ)‖ = C*(σ) e^{−2στ}` for every `τ > 0`, `1 ≤ C*(σ) < ∞`.
4. For `Re λ > −2σ` the resolvent is bounded by `C(σ)/(Re λ + 2σ)` on the whole half-plane.

This replaces the withdrawn W3 step (`docs/desk/w3_gearhart_pruss.md`, Correction) and gives a second
route to the W1 resolvent bound and the W2 uniqueness statement. Gearhart–Prüss is not used.

## 1. The identity and the conjugation

Conventions: `docs/desk/convective_weight_spectrum.md` §0. `ĝ(s) = ∫ g e^{−sz}dz`; `L²_σ` ⟺ `ĝ ∈ L²`
on `Re s = σ`; `K = e^{−z/2}T_c`, `c(s) = −tan(πs)/(2s+1)`, so `(Kg)^(s) = c(s+½)ĝ(s+½)` on
`Re s = σ−½`; `M_σ` is the maximal realization `D = {g ∈ L²_σ : −2g' + Kg ∈ L²_σ}`.

**Identity.** `N(s+½)/N(s) = [cos(πs+π/2)/cos(πs)]·Γ(2s+3)/Γ(2s+2) = −tan(πs)(2s+2)`, and
`c(s+½) = cot(πs)/(2s+2)`. The product is `−1`.

**Symbol.** `N` is analytic for `Re s > −1`; its only zeros in `|Re s| ≤ ½` are `s = ±½`. On each line
`Re s = σ'`, `|σ'| < ½`, `|N(σ'+it)| ≍ (1+|t|)^{2σ'+3/2}` (the `e^{π|t|}` of `cos` cancels the
`e^{−π|t|}` of `Γ`). The two-sided constant degenerates as `σ' → ±½` (see section 4).

**Conjugation.** Write `g = 𝒩k`, `u = e^{−z/2}`. On `Re s = σ−½`:
`(Kg)^(s) = c(s+½)N(s+½)k̂(s+½) = −N(s)k̂(s+½) = −N(s)(uk)^(s)`. So `K = 𝒩∘(−u)∘𝒩⁻¹`.
The lines `σ` and `σ−½` both lie in the open strip `(−½,½)` and the closed strip between them contains
no zero or pole of `N`, so `𝒩` and `𝒩⁻¹` are single convolution operators on both lines; the poles of
`1/N` at `±½` are never crossed. Therefore `−2g' + Kg = 𝒩(−2k' − uk)` as distributions, and
`D(M_σ) = 𝒩 D(A)`, `A = −2∂_z − e^{−z/2}`, `D(A) = {k ∈ H^p_σ : −2k' − uk ∈ H^p_σ}`.

## 2. The semigroup

`∂_τk = −2∂_zk − e^{−z/2}k` is solved along characteristics:

```
(T(τ)k)(z) = k(z−2τ) · exp(−(1−e^{−τ}) e^{−(z−2τ)/2}).
```

- Semigroup law: `(1−e^{−s}) + e^{−s}(1−e^{−t}) = 1−e^{−s−t}`.
- Bound on `H^p_σ = e^{σz}H^p(ℝ)`: translation by `2τ` has norm `e^{−2στ}`. The factor is
  `Φ(z − 2τ − 2ln x)`, `x = 1−e^{−τ}`, `Φ(z) = exp(−e^{−z/2})`. Every derivative of `Φ` is a polynomial
  in `u` times `e^{−u}`, so `Φ ∈ C_b^∞` and multiplies `H^0` and `H^3`, hence `H^p` for `p ∈ (3/2,5/2)`
  by interpolation; its multiplier norm is translation-invariant, so independent of `τ`.
- Strong continuity at `τ = 0`: `Φ(·+L) → 1` in `C^∞` on compacts as `L → ∞`, plus the uniform bound and
  density of `C_c^∞`.
- Generator: the generator `B` satisfies `B ⊆ A`. For `Re λ > −2σ`, `λ − A` is injective (the only
  distributional solution of `(λ−A)k = 0` is `e^{−λz/2}exp(+e^{−z/2})`, in no `H^p_σ`) and `λ − B` is
  onto, so `B = A`.

`S(τ) = 𝒩T(τ)𝒩⁻¹`. Because translations commute with `𝒩` and are `e^{−σh}` times a unitary, and the
damping factors at different `τ` are translates of one another, `‖S(τ)‖e^{2στ}` is the same number
`C*(σ) = ‖𝒩 e^{−u} 𝒩⁻¹‖_{L²_σ}` for every `τ > 0`.

**Closed form for a test family** (used by the probe): `k₀ = u^a e^{−bu}` gives
`Ĝ(s,τ) = N(s) e^{−2sτ} 2Γ(2s+a)(b + 1 − e^{−τ})^{−2s−a}`.

## 3. What the conjugation says about earlier statements

| Earlier statement | Status after this note |
|---|---|
| Domain condition `Ĝ(½) = 0` (Lemma D) | Automatic on `D(M_σ)`: `k̂` is analytic on `σ < Re s < σ+½`, and `N(½) = 0`. It holds as a value of the continuation; `D(M_σ) ⊄ L²_{σ+½}` (agrees with review item L1). The window `σ < ½` is this zero of `N`. |
| OBS-09 eigenmodes `G_λ` of the unconstrained operator | `k_λ = 𝒩⁻¹G_λ` has `k̂_λ = 2π²Γ(2s+λ)/sin(2πs)`, lies in `H^p_σ`, and satisfies `(λ−A)k_λ = πΓ(1+λ)·𝟙`, i.e. `(λ−M)G_λ = Γ(1+λ)·𝟙` — Lemma D's rank-one defect. Reviewer verified the constant by quadrature (2.81948734633 at `λ = 0.3`). `G_λ` is not the homogeneous solution; the implementer's first reading said it was, and that was wrong. |
| "σ(M\|L²_σ) = {Re λ = −2σ}" (noethersolve `docs/LEAP_QUEUE.md` §3) | Wrong for the maximal realization. The spectrum is the closed half-plane `Re λ ≤ −2σ`: for `Re λ < −2σ`, `E = e^{λz/2}exp(−e^{−z/2})` is a bounded functional annihilating the range of `λ − A`. Residual spectrum, no eigenvalues. Derivation only, no numerical test left of the line. |
| Lemma C: "`Re c(σ+¼) > 0` for every `σ ∈ (0,½)`" | Reviewer: false for `σ < ¼` (`sup_t Re c(σ+¼+it) < 0` at `σ = 0.01, 0.1, 0.2, 0.24`; positive from `0.26`). So the energy no-go is established only for `σ > ¼`. Not re-checked by the implementer; the form identity on the full domain was not checked by anyone. |
| Review §2: in the Theta-balance gauge the solution "leaves `L²_σ` for σ > 0" | Contradicted for `σ < ½` at the linear level: data with `ĝ(½) ≠ 0` stay in `L²_σ`, with Hölder continuity `‖S(τ)g − g‖ ~ τ^{1−2σ}` (measured exponents 0.399–0.400 at `σ = 0.3`, 0.099–0.100 at `σ = 0.45`). |
| Review item E1 (weighted transient grows with edge depth at γ = 0.5, 0.75) | Not a contradiction: those weights are outside `(0,½)`. The instant jump of `‖S(τ)‖` to `C* > 1` at `τ = 0⁺` (next section) matches E1's peak at `τ ≈ 0.003–0.05`. |
| W1 half-plane form (withdrawn) | Recovered: `‖(λ−M_σ)⁻¹‖ ≤ C(σ)/(Re λ+2σ)` for all `Re λ > −2σ`. |

## 4. Reviewer's numbers

**Resolvent against `R_slit`** (split-Cauchy formula re-implemented from
`edge_slitplane_representation.py`), worst relative difference at three strip-interior points:

| σ | λ | worst rel diff |
|---|---|---|
| 0.4 | 0.3 | 1.9e-19 |
| 0.4 | −0.5 | 2.1e-16 |
| 0.3 | −0.35+1.7i | 6.7e-15 |
| 0.2 | −0.3−4i | 1.2e-9 |
| 0.45 | −0.85+0.4i | 1.6e-7 |
| 0.3 | 2+10i | 1.5e-29 |

The 1e-7 and 1e-9 rows have `Re λ + 2σ = 0.05, 0.1`; the reviewer attributes them to an endpoint
singularity in the quadrature and did not re-run them regularised. The closed form
`R(λ)X(·;μ) = [X(·;μ) − Γ(1+μ)/Γ(1+λ)X(·;λ)]/(λ−μ)` agrees to 1e-11…1e-25 at four points with
`Re λ, Re μ > −2σ`; it does not extend left of the line.

**Lower bounds on `C*(σ)`** (sup over spans of `u^a e^{−bu}`, `ln|b| ∈ [−L, L]`, by SVD; conjugation-free
in the sense that each member's evolution was checked against the ladder equation):

| σ | L=4 | L=8 | L=16 | L=32 | converged |
|---|---|---|---|---|---|
| 0.20 | 0.999762 | 0.999996 | 1.000000 | 1.000000 | yes |
| 0.26 | 0.999829 | 0.999999 | 1.000000 | 1.000000 | yes |
| 0.28 | 1.000147 | 1.000188 | 1.000188 | 1.000188 | yes |
| 0.30 | 1.0018 | 1.0018 | 1.0018 | 1.0018 | yes |
| 0.40 | 1.0899 | 1.1066 | 1.1109 | 1.1111 | yes |
| 0.45 | 1.206167 | 1.303887 | 1.375735 | 1.394114 | no |
| 0.49 | 1.275237 | 1.484527 | 1.820081 | 2.253135 | no |
| 0.499 | 1.277877 | 1.494994 | 1.861223 | 2.413303 | no |

Readings: `C* = 1` for `σ ≤ ¼` within the families tried (contraction after removing the rate);
`C* > 1` for `σ > ¼`, so there the semigroup is not quasi-contractive (no bound `e^{ωτ}` for any `ω`);
`C*` grows as `σ → ½`. These are lower bounds. No upper bound on `C*` was computed. The reviewer's
heuristic `C* ~ (½−σ)^{−1/2}` is not fitted.

**Symbol ratio** `max/min` of `|N|/(1+|t|)^{2σ'+3/2}` over `|t| ≤ 1e4`: 4.53, 6.94, 12.2, 23.3, 113 at
`σ' = 0.10, 0.30, 0.40, 0.45, 0.49`.

**Discrete operator.** The full weighted matrix norm of `expm(τM_GG)` times `e^{2γτ}` is 5–14 at
`τ = 1` and rising; its top singular vector is a grid sawtooth near `z_min` (same phenomenon as review
item E2). On smooth Gaussian subspaces the discrete values sit below the continuum ones (γ = 0.4:
1.049–1.099 against 1.111). Two exceedances were not resolved to a mechanism: width 0.25 at `τ = 1`,
γ = 0.3 (1.119 at N = 1400, 1.055 at N = 2000); and `τ = 0.02` with the basis reaching `z_min + 1.5`
(1.16–1.19 at `z_min = −14`, absent at `−18`). The reviewer reads both as discretization or truncation;
that is a hypothesis.

## 5. Not checked by anyone

- An upper bound on `C*(σ)`; the `H^p` multiplier norm of `Φ`.
- Equality of the conjugated resolvent with `R_slit` on all of `L²_σ` (agreement is at eight parameter
  points on two data classes).
- That the manuscript's `M_σ` is literally the distributional maximal realization used here.
- Any numerical test at `Re λ < −2σ`.
- The form identity behind "dissipative for `σ < ¼`".
- The reviewer's structural remark that the kernel of `𝒩` is `−(1/2π)e^{−z}cos(e^{−z/2})`.
- The cause of the discrete sawtooth growth.
- The manuscript `.tex`: nothing in it was edited or re-read for this note.

## 6. What this does and does not move

Moves: the linear strict gap for the constrained operator (the essential half of LEM-06) from
"W1, W2 on paper, W3 withdrawn" to "one short conjugation argument, independently attacked once".
A ledger entry and a Lean statement of the symbol identity are natural next items.

Does not move: review item N1 and section 2 of the adversarial review. The nonlinear Theorem 2 lives in
Gauge (a), whose operator acts through `c` on `Re s ∈ (½, 3/2)`; conjugating it would need `1/N` across
its pole at `s = ½`, which is the rank-one defect again. Whether a stability statement for the original
system across the edge exists is still open.

## 7. Reproducing

- `python3 validate/edge_w3_conjugation_check.py` — WC-1 identity (1.5e-14), WC-2 ladder equation
  (4.3e-9), WC-3 closed form against `expm(τM_GG)`: 4.8e-4 at N = 1500, 2.4e-4 at N = 3000, against
  6.7e-2 for transport alone.
- Reviewer scripts: `validate/review_2026_11_21_w3/` (`a1` symbol; `a4*` operator-norm probe;
  `a5*` resolvent, closed form, defect constant; `sc` strong continuity; `disc*` discrete check).
  `disc.py` and `disc2.py` carry an absolute path to `src/`; `disc2.py` takes about four minutes.
