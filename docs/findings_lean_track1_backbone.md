# Lean Track 1 — the algebraic backbone of Conjecture 2.4 (DONE)

Date 2026-06-13. Formalized and machine-verified the *closed* algebraic skeleton of
the rigorous stability proof for the singular Hou–Luo self-similar profile
(Conjecture 2.4, arXiv:2604.01868). This is Track 1 of the two-track Lean plan in
`NEXT_SESSION.md` PART B — the self-contained real-analysis facts that need no
Hilbert / weighted-Sobolev / semigroup machinery. The analytic coercivity estimate
(the (★) leap) is Track 2 and is deliberately not formalized (not stubbed as a
`sorry`).

## What was proved (5 lemmas, sorry-free, zero-axiom)

Project `lean/`, pin **mathlib v4.29.0** (commit 8a178386).

| Lemma | Statement | Lean decl |
|-------|-----------|-----------|
| **L1** | steady residual of eq. (2.4) vanishes on Y>1: `R_Ω = ā·Ω̄' − c_ω·Ω̄ − Θ̄' = 0`, `R_Θ = 0` (ā=2(Y-1), Ω̄=(Y-1)^{-1/2}, Θ̄=π/2, c_l=2, c_ω=−1) | `Profile.steady_residual_Ω/Θ`, `hasDerivAt_Ωbar` |
| **L2** | edge/Θ-balance: `Ū(1)=−2` (continuous across edge), `c_l+2c_ω=0`, drift `ā(1)=0`, jump coeff `(Ū(1)+c_l)Θ̄=0` | `Profile.Ubar_one`, `theta_balance`, `drift_edge_zero`, `theta_jump_coeff_zero` |
| **L3** | c_ω-absorption: `−2(Y-1)g'(Y) − g(Y) = −2(Y-1)^{-1/2}G'(log(Y-1))` for `g=(Y-1)^{-1/2}G(log(Y-1))` — the +G/−G cancellation collapsing transport+c_ω to −2∂_z | `Profile.cω_absorption`, `hasDerivAt_gOf` |
| **L4** | weighted IBP: `∫_a^b (e^{-φ}G(-2G') + φ'e^{-φ}G²) = e^{-φ(a)}G(a)² − e^{-φ(b)}G(b)²` | `Energy.weighted_ibp` |
| **L5** | indicial roots: `−2(Y-1)g_r'−g_r = −(2r+1)(Y-1)^r`; eigen ⇔ `r=−(λ+1)/2`; admissible `r≥−1/2 ⇔ λ≤0` | `Energy.indicial_reduction`, `indicial_root_iff`, `indicial_admissible` |

These are the closed pieces the Conjecture-2.4 reduction (rigor_energy_estimate_setup.md
§7) stands on: L1–L2 make the singular profile an exact (distributional) steady state;
L3 gives the clean −2∂_z generator; L4 is the transport-dissipation half of the energy
estimate (the part that "needs no further work"); L5 is the rigorous Robin/edge
condition that drops the inadmissible homogeneous branch (the finite-codimension 𝒩
structure).

## Gates passed

- **Full build:** `lakelock lake build` → 8251 jobs, exit 0 (aggregator
  `SingularBlowup.lean` wires both files; cross-file elaboration clean).
- **Vacuity lint:** `lean_vacuity_lint.py` → PASS, 0 findings (no `True`-stubs, no
  unused named hypotheses, no subject-free structures).
- **Axiom audit:** `#print axioms` on all 11 decls → only
  `[propext, Classical.choice, Quot.sound]`. No `sorryAx`.
- **No sorry/admit/native_decide** in proof terms.

## Honest scope (what this does and does NOT do)

DOES: machine-verify the algebraic/finite skeleton — the profile genuinely solves the
steady equation, the change of variables genuinely collapses to −2∂_z, the IBP is
exact, the indicial roots are exact. These are real theorems about explicit functions,
not paraphrases.

Does NOT: prove Conjecture 2.4. The headline (THR-01) stays `speculative`; the proof
reduces to LEM-06 — the (★) bulk-coercivity of the nonlocal forcing + the H¹_W coupling
closure — which is `conjectured` (numerically de-risked in the bulk, OBS-01, but no
analytic proof). Per the no-self-certification rule, formalizing the backbone does not
raise THR-01's status. Track 2 (the (★) leap) is the genuine remaining work; it is an
analytic proof, not a near-term Lean chip (mathlib lacks the line-Hilbert /
weighted-Sobolev / C₀-semigroup infrastructure).
