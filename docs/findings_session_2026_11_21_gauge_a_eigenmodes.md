# Findings — explicit eigenmodes of the unconstrained Gauge-(a) linearization

**Date:** 2026-11-21
**Source:** independent review of `paper/stability_spectra_singular_houluo.tex` (working copy of 08:46, uncommitted on `main` at the time).
**Probe:** `validate/edge_gauge_a_eigenmode_check.py` (exit 0 on EM-1..EM-5; ~15 s, NumPy/SciPy only).
**Review form:** EM-1..EM-4 and §§1–5 as first written: single reviewer (Claude). The contour and far-field correction below and EM-5: independent Claude verification of Codex findings. Not expert-reviewed. Nothing here machine-certifies the PDE proof.
**Corrected 2026-11-21 (same day):** the first version allowed any contour `0 < s0 < 1` and said `G_λ → Γ(λ)` for every `Re λ > −2σ`. Both statements are wrong for `Re λ ≤ 0`: `Γ(2s+λ)` has a pole at `s = −λ/2`. Codex reported this with the counterexample `λ = −1/4`; I reproduced it (EM-5) and corrected §1–§3. The eigenmode claim itself survives on the full range `Re λ > −2σ`.
**Does not claim:** anything about the nonlinear PDE · the frozen-exponent operator on its constrained domain · Conjecture 2.4 · that Theorem 2 of the paper is wrong.

---

## 1. Statement

Let `X(s;λ) = Γ(2s+λ) Γ(2s+2) Γ(s) Γ(1−s)` (the paper's ladder multiplier) and

```
G_λ(z) := (1/2πi) ∫_{Re s = s0} e^{sz} X(s;λ) ds,      max(0, −Re λ/2) < s0 < 1.
```

The lower limit keeps the contour to the right of the pole of `Γ(2s+λ)` at `s = −λ/2`; `s0 = σ` is
admissible. A contour to the left of that pole gives a different function.

Then, for every `λ` with `Re λ > −2σ`, `σ ∈ (0, 1/2)`:

```
λ G_λ = −2 ∂_z G_λ + ½ Q̃_op G_λ        on all of ℝ.
```

Properties of `G_λ`:

- it lies in `L²_σ(ℝ)` for every `σ ∈ (0, 1)` with `2σ + Re λ > 0`;
- it is smooth up to the edge and vanishes linearly there: `G_λ ~ r_1 v + r_2 v² + …`, `v = e^z`,
  with `r_m = (−1)^{m+1} Γ(2m+λ) Γ(2m+2)` (from the poles of `Γ(1−s)` at `s = m`);
- far field, `z → +∞`:
  - `Re λ > 0`: `G_λ → Γ(λ)`;
  - `−2σ < Re λ ≤ 0`, `λ ≠ 0`: `G_λ = A_λ e^{−λz/2} + Γ(λ) + o(1)` with
    `A_λ = ½ Γ(2−λ) Γ(−λ/2) Γ(1+λ/2)`, the residue of `X` at `s = −λ/2`; this is still in `L²_σ`
    because `−Re λ/2 < σ`;
- `H[g](1) = −(1/π) Ĝ(1/2) = −2 Γ(1+λ) ≠ 0`.

So the operator `−2∂_z + ½Q̃_op` **without** the domain condition `Ĝ(1/2) = 0` has every `λ` in
`{Re λ > −2σ}` as an eigenvalue, including `Re λ > 0` of any size.

## 2. Why it follows from the paper's own identities

The paper's Appendix A shows `(λ+2s) X(s;λ) = c(s+½) X(s+½;λ)`, i.e. `Ĝ = X` solves the homogeneous
ladder equation. `X` is holomorphic and exponentially decaying on the strip
`max(0, −Re λ/2) < Re s < 1`, so `G_λ` has Mellin transform `X` on every line there. Since
`Re λ > −2σ > −1`, that strip contains lines `Re s = s1` and `Re s = s1 + ½` together, and Lemma B.1(1) (`𝓜[½Q̃_op G](s) = c(s+½) Ĝ(s+½)`,
stated "without assuming `H[g](1) = 0`") turns the ladder identity into the eigen-equation.

In the W2 uniqueness proof, kernel elements are `P·X` with `P` constant, and the constant is killed by
`P(0) = 0`. That step uses `(2∂_z+λ)G = e^{−z/2} T_c G ∈ L²_{σ−1/2}`, which is `½Q_op G`, not `½Q̃_op G`:
the two differ by the constant `H[g](1)/2`, whose Laplace transform on `[0,∞)` has a pole at `s = 0`.
So `P(0) = 0` holds only when `H[g](1) = 0`. Theorem 1 is therefore a statement about the operator
restricted to `{Ĝ(1/2) = 0}`, which is the frozen-exponent operator `−2∂_z + ½Q_op`.

## 3. What was checked numerically (`EM-1..EM-5`)

| Check | Content | Result |
|---|---|---|
| EM-1 | `λG + (λ+4)vG′ + 2v²G″ = ½D(v)`, with `D(v) = H[g](1+v) − H[g](1)` computed by real-space principal-value quadrature (no Mellin identity on the right-hand side); `λ = 0.3, 1, 3`, `v = 0.05, 0.3, 1, 4` | relative residual ≤ 2.2e-8 |
| EM-2 | `H[g](1) = −2Γ(1+λ)`, `G(+∞) = Γ(λ)`, contour independence | agree to 1e-6 (the far-field limit is asserted for `λ ≥ 1` only; at `λ = 0.3` the `e^{−λz/2}` correction is still visible at `z = 40`) |
| EM-3 | edge coefficients `r_1`, `r_2` against the residue formula; `|r_{m+1}/r_m|` grows like `(2m)⁴` | agree to 0.2% |
| EM-4 | the repo's discretized Gauge-(a) matrix (`src/stability_edge_z.py`, as used in `edge_nonlinear_comoving_closure.py`) | no unstable eigenvalue; marching `G_1` gives `H[g](1)(T)/H[g](1)(0) ≈ e^{−T}`, where the continuum mode has `e^{+T}` |
| EM-5 | negative `λ = −1/4` (pole at `s = 1/8`): contours `s0 = 0.25` and `0.9` agree; `s0 = 0.1` differs by exactly `A_λ e^{z/8}`; far field `A_λ e^{z/8} + Γ(λ)` with remainder matching the next pole term; eigen-identity; `H[g](1)` | `A_λ = 4.6506249067`; `G(16) = 29.4658775645` (Codex's value), `G(32) = 249.01` against `Γ(λ) = −4.90`; eigen-identity residual ≤ 1.6e-8; `H[g](1)` agrees to 1e-6 |

EM-1..EM-3 use positive `λ` only (`0.3, 1, 3`). EM-5 is the only negative-`λ` test, and only real `λ`
is tested anywhere in the probe.

The identity in EM-1 is the `v`-derivative of `v·[λG + 2vG′ − ½Q̃_op G] = 0`; both sides vanish at
`v → 0`, so it is equivalent to the eigen-equation.

## 4. Consequences for the paper

1. **Theorem 1 wording.** "Linear spectral gap" / "absence of spectrum in `{Re λ > −2σ}`" is true only
   on the constrained domain. The Gauge-(a) operator that Theorem 2 evolves carries no such
   constraint, and in `L²_σ(ℝ)` its point spectrum contains the whole half-plane. The paper should
   say this next to Theorem 1.
2. **No contradiction with Theorem 2.** The Taylor series of `G_λ` at the edge has zero radius of
   convergence (EM-3), so `G_λ ∉ 𝓗¹(D(1,1))`. Edge analyticity is doing essential work.
3. **The "open" `C^k` extension.** `G_λ` is `C^∞` up to the edge with `G = O(v)`. In this gauge the
   linearized equation therefore has solutions `e^{λτ} G_λ` of every growth rate in that class. As
   posed, the extension faces explicit growing modes, not only a missing estimate.
4. **Numerical Theorem.** The discretized operator does not carry these modes (EM-4). The grid stops at
   `z_min`, so every grid function vanishes identically near the edge and the edge term in
   `d/dτ H[g](1)` (the `−c_1/3` term in `a′ = −a − c_1/3`) is absent. `#{Re λ > 0} = 0` on the grid is
   evidence about the truncated problem, not about the continuum operator on non-analytic edge data.

## 5. Open questions (not settled here)

- Whether these modes are physically relevant or an artifact of fixing the edge amplitude through the
  unbounded functional `H[g](1)` (a different gauge might remove them).
- Whether the continuum Gauge-(a) evolution of grid-like data (zero near the edge) stays close to the
  discretized evolution over the times the gates use.
- Mechanism in words, unverified beyond a leading-order estimate: for `G ≈ T` outside an edge layer of
  width `W`, `H[g](1) ≈ −(4T/π) W^{−1/2}` and the far field obeys `Ṫ = −H[g](1)/2`, a positive feedback
  with rate `~ W^{−1/2}`.
