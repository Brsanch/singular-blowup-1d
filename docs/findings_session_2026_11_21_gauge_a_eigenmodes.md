# Findings — explicit eigenmodes of the unconstrained Gauge-(a) linearization

**Date:** 2026-11-21
**Source:** independent review of `paper/stability_spectra_singular_houluo.tex` (working copy of 08:46, uncommitted on `main` at the time).
**Probe:** `validate/edge_gauge_a_eigenmode_check.py` (exit 0 on EM-1..EM-4; ~10 s, NumPy/SciPy only).
**Review form:** single reviewer (Claude). Not independently audited, not expert-reviewed.
**Does not claim:** anything about the nonlinear PDE · the frozen-exponent operator on its constrained domain · Conjecture 2.4 · that Theorem 2 of the paper is wrong.

---

## 1. Statement

Let `X(s;λ) = Γ(2s+λ) Γ(2s+2) Γ(s) Γ(1−s)` (the paper's ladder multiplier) and

```
G_λ(z) := (1/2πi) ∫_{Re s = s0} e^{sz} X(s;λ) ds,      0 < s0 < 1.
```

Then, for every `λ` with `Re λ > −2σ`, `σ ∈ (0, 1/2)`:

```
λ G_λ = −2 ∂_z G_λ + ½ Q̃_op G_λ        on all of ℝ.
```

Properties of `G_λ`:

- it lies in `L²_σ(ℝ)` for every `σ ∈ (0, 1)` with `2σ + Re λ > 0`;
- it is smooth up to the edge and vanishes linearly there: `G_λ ~ r_1 v + r_2 v² + …`, `v = e^z`,
  with `r_m = (−1)^{m+1} Γ(2m+λ) Γ(2m+2)`;
- it tends to the constant `Γ(λ)` as `z → +∞`;
- `H[g](1) = −(1/π) Ĝ(1/2) = −2 Γ(1+λ) ≠ 0`.

So the operator `−2∂_z + ½Q̃_op` **without** the domain condition `Ĝ(1/2) = 0` has every `λ` in
`{Re λ > −2σ}` as an eigenvalue, including `Re λ > 0` of any size.

## 2. Why it follows from the paper's own identities

The paper's Appendix A shows `(λ+2s) X(s;λ) = c(s+½) X(s+½;λ)`, i.e. `Ĝ = X` solves the homogeneous
ladder equation. `X` is holomorphic and exponentially decaying on the strip `0 < Re s < 1`, so `G_λ`
has Mellin transform `X` on every line there, and Lemma B.1(1) (`𝓜[½Q̃_op G](s) = c(s+½) Ĝ(s+½)`,
stated "without assuming `H[g](1) = 0`") turns the ladder identity into the eigen-equation.

In the W2 uniqueness proof, kernel elements are `P·X` with `P` constant, and the constant is killed by
`P(0) = 0`. That step uses `(2∂_z+λ)G = e^{−z/2} T_c G ∈ L²_{σ−1/2}`, which is `½Q_op G`, not `½Q̃_op G`:
the two differ by the constant `H[g](1)/2`, whose Laplace transform on `[0,∞)` has a pole at `s = 0`.
So `P(0) = 0` holds only when `H[g](1) = 0`. Theorem 1 is therefore a statement about the operator
restricted to `{Ĝ(1/2) = 0}`, which is the frozen-exponent operator `−2∂_z + ½Q_op`.

## 3. What was checked numerically (`EM-1..EM-4`)

| Check | Content | Result |
|---|---|---|
| EM-1 | `λG + (λ+4)vG′ + 2v²G″ = ½D(v)`, with `D(v) = H[g](1+v) − H[g](1)` computed by real-space principal-value quadrature (no Mellin identity on the right-hand side); `λ = 0.3, 1, 3`, `v = 0.05, 0.3, 1, 4` | relative residual ≤ 2.2e-8 |
| EM-2 | `H[g](1) = −2Γ(1+λ)`, `G(+∞) = Γ(λ)`, contour independence | agree to 1e-6 (the far-field limit is asserted for `λ ≥ 1` only; at `λ = 0.3` the `e^{−λz/2}` correction is still visible at `z = 40`) |
| EM-3 | edge coefficients `r_1`, `r_2` against the residue formula; `|r_{m+1}/r_m|` grows like `(2m)⁴` | agree to 0.2% |
| EM-4 | the repo's discretized Gauge-(a) matrix (`src/stability_edge_z.py`, as used in `edge_nonlinear_comoving_closure.py`) | no unstable eigenvalue; marching `G_1` gives `H[g](1)(T)/H[g](1)(0) ≈ e^{−T}`, where the continuum mode has `e^{+T}` |

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
