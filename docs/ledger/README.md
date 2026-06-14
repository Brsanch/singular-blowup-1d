# Proof-obligation ledger — singular-blowup-1d

Output discipline for the Conjecture-2.4 (linear stability of the singular
Hou–Luo profile) proof arc. Schema: `../../../references/proof_obligation_ledger_schema.md`.

One file per claim. Status taxonomy is no-self-certification: `verified` means
Lean-proved zero-sorry/zero-axiom **per entry** (or a `validate/` gate that
re-derives independently); the headline THR caps at `speculative` until external
review, no matter how green the lemma DAG is.

## Current map

| id | claim | status | machine check |
|----|-------|--------|---------------|
| THR-01 | Conjecture 2.4: singular profile is asymp. stable mod 𝒩 | speculative (0.70) | — |
| LEM-01 | steady residual of eq. (2.4) vanishes on Y>1 | verified | Lean + validate gate |
| LEM-02 | edge/Θ-balance: U(1)=−2, c_l+2c_ω=0, ā(1)=0 | verified | Lean |
| LEM-03 | c_ω-absorption: transport+c_ω ⇒ −2∂_z | verified | Lean (+1e-11 gate) |
| LEM-04 | weighted IBP transport-dissipation identity | verified | Lean |
| LEM-05 | edge indicial roots r=−(λ+1)/2; admissible ⇔ Reλ≤0 | verified | Lean |
| LEM-06 | (★) bulk coercivity of the nonlocal forcing | conjectured (0.60) | validate gate (bulk only) |
| OBS-01 | bulk energy abscissa < 0 (anisotropic weight) | — | validate/energy_estimate_search.py |
| OBS-02 | steady residual ~2e-16; U(1)=−2 exact | — | validate/verify_singular_steady.py |

LEM-01..05 = the **Track-1 algebraic backbone** (Lean-verified). LEM-06 (the (★)
leap) + the H¹_W coupling = **Track 2**, the genuine remaining analytic work; it
stays `conjectured` with the OBS gates as evidence. THR-01 is `speculative`.
