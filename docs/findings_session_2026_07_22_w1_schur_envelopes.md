# Findings — W1 quantitative Schur / Plemelj envelopes (2026-07-22)

**Status:** WIN (new numbers + one desk-prose correction).  
**Gate:** `validate/edge_schur_envelopes.py` (pure core `src/schur_envelopes.py`).  
**Dual-run:** bit-identical KEY line (rel tol 0 within float print).  
**Does not claim:** full W1 write-up / LEM-06 closed / semigroup.

## Claim attempted

Get explicit numerical constants for the dressed two-weight Schur kernel of desk
§10 as functions of `σ` and `δ = Re λ + 2σ` (off-diagonal envelope, effective
decay rate, diagonal-band max, termwise `C`, defect-density rate).

## Novelty vs prior art in-repo

| Prior (desk §10 / OBS-07/08) | This session |
|---|---|
| Spot checks: `dressed·e^{π\|Δτ\|}` = 0.117→0.022 etc. at fixed `(σ=0.4,λ=0.3)` | Full table of `S_off(σ,δ)=sup_{\|Δτ\|≥1} dressed·e^{π\|Δτ\|}` |
| Rate stated as `≥ e^{−π\|Δτ\|}` (margin only) | Measured `α_eff` on former-corner wing |
| Termwise `≈1.15/dist` (one bulk remark) | `C_term(σ,δ)=(Σ_n sup\|Π_n|)·δ` table |
| Defect kernel “`~e^{−3π\|τ\|/2}` decaying” (prose) | Bare line density **grows** `γ≈3π` |

No prior `S_off`, `α_eff`, `M_diag(δ)`, or `γ_def` numbers appear in `docs/findings_*` or ledger.

## Key numbers (anchor `σ=0.4`, `δ=1.1` ⇔ `λ=0.3`; tol **1e-3** rel)

```
S_off_anchor          = 0.551085
alpha_eff_anchor      = 3.71162     # π + 0.570  (former-corner wing, τ'=-12)
C_term_anchor         = 1.12545     # desk ~1.15 class, tightened
M_diag_anchor         = 10.2032     # max dressed on |Δτ|≤1, τ≠τ'
defect_growth_anchor  = 8.89354     # ≈ 3π (=9.4248); GROWTH not decay
S_off × δ (δ=0.15)    = 0.12414
```

Gate KEY line (both runs identical):

```
KEY S_off_anchor=0.551085 alpha_anchor=3.71162 C_term_anchor=1.12545
    M_diag_anchor=10.2032 defect_growth_anchor=8.89354 S_off_x_delta_small=0.12414
```

## Full table (band=1, `n_τ=41`, dps=20)

| σ | δ | λ | S_off | α_eff | M_diag | C_term | γ_def |
|---|---|---|-------|-------|--------|--------|-------|
| 0.25 | 0.15 | −0.350 | 1.049 | 3.4647 | 28.04 | 1.3975 | 9.0630 |
| 0.25 | 0.30 | −0.200 | 0.8475 | 3.4949 | 23.11 | 1.3486 | 9.0425 |
| 0.25 | 0.50 | 0.000 | 0.7478 | 3.5351 | 17.43 | 1.2992 | 9.0151 |
| 0.25 | 1.10 | 0.600 | 0.6029 | 3.6558 | 9.276 | 1.2090 | 8.9336 |
| 0.25 | 1.50 | 1.000 | 0.5297 | 3.7354 | 6.983 | 1.1737 | 8.8797 |
| 0.40 | 0.15 | −0.650 | 0.8276 | 3.5205 | 30.84 | 1.1398 | 9.0230 |
| 0.40 | 0.30 | −0.500 | 0.7294 | 3.5507 | 25.42 | 1.1335 | 9.0024 |
| 0.40 | 0.50 | −0.300 | 0.6795 | 3.5910 | 19.17 | 1.1358 | 8.9751 |
| 0.40 | 1.10 | 0.300 | 0.5511 | 3.7116 | 10.20 | 1.1254 | 8.8935 |
| 0.40 | 1.50 | 0.700 | 0.5048 | 3.7912 | 7.681 | 1.1117 | 8.8397 |

## Interpretation (firewall-clean)

1. **Off-diagonal Schur is quantitatively closed on the grid.**
   `S_off ≤ 1.05` over the probed `(σ,δ)`, and `S_off·δ ∈ [0.12, 0.80]`
   (product ratio ~6 across a 10× δ-range). The π-margin is nearly sharp on
   the former corner: `α_eff − π ∈ [0.32, 0.65]`, not a large slack.

2. **Termwise bulk constant is δ-mild.**
   `C_term ∈ [1.11, 1.40]`; at the desk anchor `1.125` (was “≈1.15”).
   Supports `‖R‖ ≲ C/δ` with `C ≲ 1.4` on this termwise proxy (not yet the
   full operator norm).

3. **Desk §10 defect-kernel prose needs correction (new obstruction to a
   naive W1 write-up step).**
   Bare density `k(τ)=1/(|X(σ+iτ)| |λ+2s|)` **grows** as `e^{γ|τ|}` with
   `γ ∈ [8.84, 9.06] ≈ 3π` (Stirling: `|X|~e^{−3π|τ|}`). It does **not**
   decay as `e^{−3π|τ|/2}`. W1 must either (i) identify a different functional
   kernel (w-plane Jacobian / split measure / restricted Fhat class) that is
   L²-bounded, or (ii) bound the defect via the ladder sum rather than a
   decaying line density. This is a **named bookkeeping correction**, not a
   refutation of the slit-plane representation (P1 of `edge_slitplane_representation.py`
   still holds).

## Gate PASS criteria (all green, dual-run)

- P1 `α_eff ≥ π − 0.05` all rows  
- P2 `S_off` finite, `S_off·δ < 50`  
- P3 `C_term ∈ [0.5, 3]` at desk anchor  
- P4 `S_off` grows as `δ↓`; product ratio < 50  
- P5 `γ_def ∈ (2.5π, 3.5π)` (growth, not decay)

## What this does *not* close

- Full W1 rigorous write-up (Plemelj BVP + Stirling bookkeeping with the
  corrected defect kernel).  
- Lemma E′ / W2, Gearhart–Prüss / W3.  
- `M_Gh` coupling and nonlinear legs.  
- Operator-norm `‖R(λ)‖` (only termwise + Schur-envelope proxies).

## Reproduce

```bash
python3 validate/edge_schur_envelopes.py
# KEY line must match within 1e-3 rel
```
