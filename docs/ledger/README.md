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
| LEM-06 | (★) coercivity — reduced to a Wiener–Hopf index (LEM-07) | conjectured (0.60) | validate gates (bulk only) |
| LEM-07 | edge Mellin symbol 𝓗(p)=cot(πp); forcing −½ ladder; ⟹ Mellin difference eqn | proved-on-paper | validate/edge_mellin_symbol.py + mpmath/sympy (4 methods) |
| OBS-01 | bulk energy abscissa < 0 (anisotropic weight) | — | validate/energy_estimate_search.py |
| OBS-02 | steady residual ~2e-16; U(1)=−2 exact | — | validate/verify_singular_steady.py |
| OBS-03 | +40 edge abscissa = inadmissible-branch resonance (energy route dead) | — | validate/edge_layer_obstruction.py |
| OBS-04 | no Re>0 mode, but max Re = marginal essential spectrum on Re=0 (not a gap) | — | validate/edge_spectrum_continuation.py |
| OBS-05 | arg-principle count = 0 (in-strip corroboration; NOT the WH index) | — | validate/edge_argprinciple_count.py |
| OBS-06 | Mellin-line: lattice solns edge-divergent ∀λ; winding = indicial tautology (NOT the index) | — | validate/edge_wh_structure.py |

LEM-01..05 = the **Track-1 algebraic backbone** (Lean-verified). LEM-06 (the (★)
leap) + the H¹_W coupling = **Track 2**, the genuine remaining analytic work; it
stays `conjectured` with the OBS gates as evidence. THR-01 is `speculative`.

**Edge-layer attack (2026-06-13, `../desk/edge_layer_estimate.md`):** 6-lens derive +
18 adversarial verifies + independent numerical reproduction (OBS-03) → the edge leg of
LEM-06 does NOT close by any weighted-energy/Schur/flux method (proven dead, not just
undone). The leap is rescoped to a non-self-adjoint **spectral/Fredholm** bound (resolvent
or Mellin/Wiener–Hopf symbol at `r=−(λ+1)/2`). Circularity firewall: `H(g)(1)=0` is
automatic for L²_W modes — must NOT be asserted from the benign spectrum (NS-C6/SQG trap).
