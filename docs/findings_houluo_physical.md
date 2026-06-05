# Finding: physical-DNS confirms c_ω=−1, but c_l needs dynamic rescaling

Date: 2026-06-05.  Solver: `src/houluo.py` (physical-space pseudospectral,
box [−π,π), 2/3-dealias, RK4).  Model `ω_t+uω_x=θ_x`, `θ_t+uθ_x=0`, `u_x=H[ω]`.

## Result

Finite-time focusing blowup confirmed for smooth (non-degenerate) ICs, focusing
at a symmetry point (x≈0 or x≈±π).  The amplitude exponent is recovered cleanly
and grid-robustly:

| IC (ω₀, θ₀) | focus x* | fitted c_ω (expect −1) |
|---|---|---|
| (sin x, −cos x) | ≈ π | **−1.0033** |
| (sin x, 1−cos x) | ≈ π | **−0.9959** |
| (sin x, cos x) | ≈ 0 | −1.157 (messy multi-focus transient) |

`max|ω| ~ (T−t)^{−1}` to <0.4% on the clean cases — matches 2308.01528 (c_ω=−1).

## Why c_l is NOT obtainable from physical DNS (expected, not a bug)

The fitted `c_l` came out ≈ −0.03…−0.07, nowhere near the predicted (2, 4.53).
Reason: the length scale `L(t)=max|ω|/max|ω_x|` collapses and then **stalls at
the grid scale** (L ≈ 0.0022 ≈ 0.75·dx at N=2048; ≈ 1.1e-3 ≈ 0.72·dx at N=4096).

Because `c_ω=−1` while `c_l∈(2,4.53)`, the self-similar core width `(T−t)^{c_l}`
shrinks **far faster** than the amplitude `(T−t)^{−1}` grows.  Per amplitude
doubling, (T−t) halves and the core shrinks by `2^{c_l}` = 4…23×.  So the core
goes **sub-grid almost immediately**, and doubling N (2048→4096) barely extends
the resolved window — confirmed: the resolution alarm fires *earlier* (at lower
|ω|) at N=4096.  This is the same sub-grid-scale lesson as the NS regularity
work: a scale `~(T−t)^{c_l}` with `c_l>1` is unreachable at fixed grid.

**Consequence:** c_l, the profile shape, the far-field decay, and the
smooth-vs-singular (degenerate-IC) distinction all require the
**dynamic-rescaling formulation** (eq 2.4 of 2604.01868), which holds the
profile at O(1) scale on the grid for all rescaled time.  That is the next build
and the real workhorse; its linearization is also the stability operator
(Conjecture 2.4 / task 7).

## Diagnostic bug fixed

`tail_fraction` originally measured |k|>N/3, which the 2/3 dealias mask zeroes
every step → it read ~1e-26 (useless) even as the core hit grid scale.  Fixed to
the upper RETAINED band |k|∈[N/6,N/3]; it now correctly fires at L≈grid scale.
