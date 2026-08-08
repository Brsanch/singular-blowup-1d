# Findings — W1-D1 reopen (cancellation outlines + dual probe)

**Date:** 2026-08-07  
**Status:** OUTLINES + NUMERIC-FORECAST only. **Does not claim D1 / LEM-06.**  
**Artifacts:**
- `docs/desk/w1_d1_cancellation_outlines.md`
- `validate/edge_d1_dual_probe.py`

## Why this session

Magnet #4 (2026-07-25) parked on leap-shaped wall **W1-D1**. Reopen rule: only with a
cancellation-lemma *candidate*, not more Schur/OBS-07 grinding. This session wrote
ranked outlines and a dual probe.

## Probe results (lam=0.3, σ=0.4, δ=1.1 unless noted)

### P0 — ladder coefficients
| k | \|a_k\| | sign |
|---|---|---|
| 0 | 1.36e-1 | + |
| 1 | 1.50e-3 | − |
| 2 | 2.19e-6 | + |
| … | super-factorial decay | (−1)^k |

`Σ a_k = C_∞ ≈ 0.1349 ≠ 0`.

### P1 — empirical R = |D1(F)| / ‖F‖_{L²(line)} on bilateral bumps
`F(s)=2 K_s(1) e^{-s z_0}`:

| z₀ | R |
|---|---|
| +2 | 0.201 |
| 0 | 0.241 |
| −2 | 0.254 |
| −4 | 0.037 |
| −6 | 0.136 |
| −8 | **1.63** |

**Correction to a naive OBS-07 reading:** edge-localized data can make **|D1| large**
(R ≳ 1) even though the *full* resolvent gain ‖G‖/‖F‖ annihilates with depth
(OBS-07). Mechanism: at depth `C₀(s)≈C₀(½)`, so `G=X[C₀−C₀(½)]` is small while the
scalar `C₀(½)=D1` need not be. Two-stage cancellation is load-bearing; bounding D1
is **not** implied by small gains.

### P2 — formal H² dual ψ
```
ψ(τ) = (1/2π) Σ_k a_k / (σ + iτ − s_k),  s_k = ½+k
```
Truncated ‖ψ‖_{L²([-T,T])} **saturates** (~0.121 at T=5…200). Triangle bound
`Σ |a_k| ‖K_k‖₂/(2π) ≈ 0.122`.  

**Correction:** earlier outline draft claimed `1/τ ∉ L²` at infinity — **false**
(`∫_{|τ|>1} τ^{-2}dτ < ∞`; singularity of `1/τ` is at 0, where ψ is regular).
So **ψ ∈ L²(ℝ)** as an abstract H² representer.

### P3 — δ-scaling of R at z₀=0
| δ | R | R·δ |
|---|---|---|
| 1.5 | 0.183 | 0.275 |
| 1.1 | 0.241 | 0.265 |
| 0.6 | 0.299 | 0.179 |
| 0.3 | 0.311 | 0.093 |
| 0.15 | 0.306 | 0.046 |

On this family **R is O(1) in δ**, not O(1/δ). Either (i) the wall’s `O(1/δ)` target
is for the full resolvent not the scalar D1, or (ii) bumps do not saturate the
worst-case D1, or (iii) the target should be revised. **Do not cash O(1/δ) for D1
from this family.**

### Cauchy recovery kill of naive Outline A on L²_σ data
For the same bilateral bumps, the H² right-half Cauchy formula
```
(1/2π) ∫ F(σ+iτ)/(σ+iτ − s_k) dτ
```
**does not recover** `F(s_k)` (relative errors O(1)–O(2) at k=0,1,2; z₀=0,−4,+2).
Bilateral `e^{-cosh(z-z₀)}` is **not** in `H²(Re s > σ)`. Plancherel for `L²_σ`
is full `L²(line)`, not Hardy.

**Conclusion:** Outline A gives a continuous dual on the **Hardy subspace only**.
It does **not** define or bound D1 on the physical `L²_σ` line data used by the
ladder/resolvent construction. Not a discharge path without a new identification
of the relevant subspace (or a different integral formula).

## Outline status after probe

| ID | Route | Status |
|---|---|---|
| A | H²-Cauchy dual ψ | ψ∈L² OK; **Cauchy recovery fails on bilateral dense class** → not D1 for L²_σ |
| B | Coboundary telescope | Reformulation; still needs multi-line bounds |
| C | Bent Mellin–Barnes | Still the named classical full-resolvent route (desk §8) |
| D | Slit-plane kernel for D(1) | Best next dual probe: large-τ expansion of the *scalar* slit kernel |

## What not to do next
- Re-run Schur envelopes / OBS-07 as if they bound D1.
- Claim D1 from P1 ratios (family lower bounds only; edge R already >1).
- Open CAP on D1 (no estimate outline worth certifying yet).

## Honest next chips (pick one)
1. **Outline D probe:** extract the integral kernel of slit-plane `D(1)` against
   line `F` (from `edge_slitplane_representation.py` identities) and test
   `K_slit ∈ L²` numerically on a τ-grid; compare large-τ coeff to `C_∞`.
2. **Outline C execution slice:** bent-contour Barnes for the *scalar* `C₀(½)`
   only (smaller than full resolvent Schur).
3. **Ship path:** freeze linear rigor; arXiv the ready numerical spectra paper
   (`paper/stability_spectra_singular_houluo.*`) — magnet #4 option 1.

## What I did NOT verify
- Full operator norm of D1 (only bump-family Rayleigh quotients).
- Whether a left/right split of bilateral data recovers F via two Hardy projections.
- Muskhelishvili uniqueness (R1) or continuum E′ / GP.
- Paper abstract metadata / arXiv category fit for upload.
