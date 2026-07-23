# Findings — Lemma E′ residual floor (2026-07-22)

**Status:** NEW QUANTITATIVE EVIDENCE for E′ (not a proof).  
**Gate:** `validate/edge_eprime_probe.py`  
**Dual-run:** KEY bit-identical (rel 0 within print).

## Claim

Is an L²_σ “eigenfunction” without strip-regularity constructible from the
homogeneous edge branch `G ~ e^{−λz/2}`, or does residual force the bootstrap?

## Result (new numbers)

At `σ=0.4`, trial `λ=0.3` (above transport line `−0.8`):

| Quantity | Value | Meaning |
|---|---|---|
| **γ_w** | **1.10035** | edge-mass growth rate of `G_hom` in L²_σ (= `Reλ+2σ` to 0.03%) |
| **r_min** | **0.58576** | min over depths of `‖(λ−M)G_hom‖_σ/‖G_hom‖_σ` |
| **r(zmin=−14)** | **32.583** | residual **grows** with depth (ratio r(−14)/r(−6)=**55.6**) |
| **λ_spur (free edge)** | **+1.7066** | unstable mode when `edge_bc=False` |
| **α_edge free** | **−0.544** | homogeneous-class growth (not recessive) |
| **max Re (edge_bc)** | **−0.00156** | no unstable mode with inflow BC |
| **α_edge BC** | **+3.03** | rightmost BC mode not homogeneous-growing |

KEY line (both runs):
```
KEY gamma_w=1.10035 r_min=0.585763 lam_spur_re=1.7066 alpha_free=-0.544389 alpha_bc=3.02852 maxRe_bc=-0.00156221
```

## Interpretation

1. **No constructible approximate eigenfunction from `G_hom`.**  
   Residual floor `r ≥ 0.59` on the whole depth grid; deeper truncations make residual
   *worse* (nonlocal `K` acting on edge mass), not better. A counterexample to E′ would
   need `r → 0` as `zmin → −∞`. Observed: opposite.

2. **Weight divergence is exact.**  
   `γ_w = Reλ+2σ` to 0.03% — the homogeneous branch is not in L²_σ for `Reλ > −2σ`,
   measured not just argued.

3. **Free-edge unstable mode is the inadmissible branch.**  
   `Re λ ≈ +1.71` with edge growth `α < 0` (homogeneous-class). Imposing
   `G(z_min)=0` kills it and restores `max Re ≤ 0`. This pins the discrete inflow BC
   as the same object E′ names in the continuum.

## What this is *not*

- Not a proof of strip-regularity / E′.  
- Not a continuum Mellin argument (discrete operator residual).  
- Does not replace the analytic Liouville step of Prop. E.

## Reproduce

```bash
python3 validate/edge_eprime_probe.py
```
