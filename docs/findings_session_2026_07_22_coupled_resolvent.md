# Findings — Coupled (G,h) resolvent scan (2026-07-22)

**Status:** NEW GAIN TABLE for Θ-drive / coupled forcing (not a pure-Ω rehash of OBS-07).  
**Gate:** `validate/edge_coupled_resolvent_scan.py`  
**Dual-run:** KEY bit-identical.

## Question

OBS-07: pure-Ω continuum ladder has edge-forcing gain `∝ e^{z₀/2}` (operator annihilates
edge input). Open: does the `M_Gh = e^{−z/2} ∂_z` coupling reintroduce edge growth when
Θ is forced at the edge? Does the full block-triangular operator keep the gap?

## Setup

Frozen discrete edge-resolved operator (`src/stability_edge_z.py`), `N=900`,
`z∈[−14,8]`, `edge_bc=True`, `σ=0.4`, `λ=0.3`.

Block-triangular solve:
```
(λ − M_hh) h = f
(λ − M_GG) G = F + M_Gh h
```

Forcings: bump at `z₀ ∈ {−2,−4,−6,−8,−10}` —
(A) pure Ω (`f=0`), (B) pure Θ (`F=0`), (C) equal coupled.

## Gain table (weighted L²_σ)

| z₀ | gain_Ω | gain_Θ→G | gain_coupled | ‖h‖_Θ/‖f‖ |
|----|--------|----------|--------------|-----------|
| −2 | 0.3035 | 0.4298 | 0.6923 | 0.663 |
| −4 | 0.1721 | 0.6522 | 0.7402 | 0.663 |
| −6 | 0.0821 | 0.8382 | 0.7962 | 0.663 |
| −8 | 0.0345 | 0.9479 | 0.8349 | 0.663 |
| −10 | 0.0133 | **0.9892** | 0.8486 | 0.663 |

Log-slopes `d log(gain)/dz₀`:

| channel | slope | meaning |
|---------|-------|---------|
| pure Ω | **+0.393** | annihilation class (OBS-07: → +½) |
| Θ→G | **−0.102** | mild growth, saturates to O(1) |
| coupled | **−0.026** | essentially flat O(1) |

KEY (both runs):
```
KEY maxRe=-1.66802e-06 slopeO=0.392767 slopeThG=-0.10206
    gO_m2=0.30349 gO_m10=0.013348 gThG_m2=0.429774 gThG_m10=0.98919 ratio_deep=74.1078
```

## New conclusions

1. **Coupling does not reintroduce Re≥0 spectrum.**  
   `max Re(full coupled) = −1.7e−6`, `#Re>0 = 0`.

2. **`M_Gh` does not destroy edge control.**  
   Despite `M_Gh ~ e^{−z/2}` blowing at the edge, the Θ→G gain **saturates to ≈0.99**
   (growth factor only 2.3 from z₀=−2 to −10). No `e^{−z₀/2}` blowup of the response.

3. **Channel asymmetry is structural.**  
   Pure-Ω edge forcing is annihilated (gain ÷23 from −2 to −10).  
   Pure-Θ edge forcing drives an O(1) G-response.  
   Deep ratio `(Θ→G)/Ω ≈ 74` is **Ω-annihilation**, not Θ blowup — do not read it as
   a coupling instability.

4. **Θ-sector resolvent is depth-independent.**  
   `‖h‖_σ/‖f‖_σ = 0.663` flat in z₀ (pure transport `M_hh = −2∂_z`).

## Novelty vs prior art

| Prior | This |
|-------|------|
| OBS-07 pure-Ω continuum ladder gains | Θ→G and coupled discrete gains |
| `theta_munu_gap.py` spectrum/decay rates | resolvent gain-vs-depth table |
| coevolve qualitative “no coupling runaway” | quantitative slope_Θ→G = −0.102, gThG→0.99 |

## What this is *not*

- Not continuum Mellin/ladder for the coupled system (discrete Hilbert remains).  
- Not a proof that `‖R_coupled‖ ≲ C/dist`.  
- Not nonlinear / co-moving.

## Reproduce

```bash
python3 validate/edge_coupled_resolvent_scan.py
```
