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
| LEM-06 | (★) — point half CLOSED (LEM-08); essential half: PROOF SKELETON COMPLETE via the slit-plane Plemelj representation 𝒢 = X·[D(w(s))−D(1)] + dressed-Schur (π-margin, all quadrants); remaining = W1 write-up, W2 Lemma E′, W3 generation+GP | conjectured (0.75) | LEM-09 gate + OBS-07/08 + edge_slitplane_representation.py |
| LEM-07 | edge Mellin symbol 𝓗(p)=cot(πp); forcing −½ ladder; ⟹ Mellin difference eqn | proved-on-paper | validate/edge_mellin_symbol.py + mpmath/sympy (4 methods) |
| LEM-08 | no Re λ≥0 POINT eigenvalue mod 𝒩 (exact integer obstruction; τ elementary) | proved-on-paper | validate/edge_pointspectrum_rigid.py (+ pole_collision, route3, firewall_audit) |
| LEM-09 | convective-weight structure pack: ε-similarity; transport line ⊆ σ_ess; energy no-go w/ exact σ+¼-band mechanism; rank-one defect = domain condition; closed-form defect-closed ladder resolvent (index-0) | proved-on-paper (E-part mod Lemma E′) | validate/edge_convective_structure.py (12 checks) |
| OBS-01 | bulk energy abscissa < 0 (anisotropic weight) | — | validate/energy_estimate_search.py |
| OBS-02 | steady residual ~2e-16; U(1)=−2 exact | — | validate/verify_singular_steady.py |
| OBS-03 | +40 edge abscissa = inadmissible-branch resonance (energy route dead) | — | validate/edge_layer_obstruction.py |
| OBS-04 | no Re>0 mode, but max Re = marginal essential spectrum on Re=0 (not a gap) | — | validate/edge_spectrum_continuation.py |
| OBS-05 | arg-principle count = 0 (in-strip corroboration; NOT the WH index) | — | validate/edge_argprinciple_count.py |
| OBS-06 | Mellin-line: lattice solns edge-divergent ∀λ; winding = indicial tautology (NOT the index) | — | validate/edge_wh_structure.py |
| OBS-07 | edge-forcing resolvent gain ∝ e^{z₀/2} = ε, λ-uniform (edge input annihilated; supports σ(M|_σ) = transport line, ‖R‖~C/dist) | measured (0.85) | validate/edge_convective_resolvent_scan.py |
| OBS-08 | high-frequency corner benign: dressing → 1 with |Ω| at every depth (the ex-corner of the crude Schur; coupling c(σ+iΩ)~1/(2|Ω|)) | measured (0.85) | validate/edge_slitplane_representation.py (P3) |
| OBS-09 | Gauge-(a) operator WITHOUT the domain condition Ĝ(½)=0 has an explicit eigenfunction G_λ = 𝓜⁻¹[X(·;λ)] for every Re λ > −2σ (C^∞ at the edge, O(v), not edge-analytic, H[g](1) = −2Γ(1+λ) ≠ 0); the truncated grid does not carry it | measured (0.9) | validate/edge_gauge_a_eigenmode_check.py (EM-1..EM-5) |

LEM-01..05 = the **Track-1 algebraic backbone** (Lean-verified). LEM-06 (the (★)
leap) + the H¹_W coupling = **Track 2**, the genuine remaining analytic work; it
stays `conjectured` with the OBS gates as evidence. THR-01 is `speculative`.

**Edge-layer attack (2026-06-13, `../desk/edge_layer_estimate.md`):** 6-lens derive +
18 adversarial verifies + independent numerical reproduction (OBS-03) → the edge leg of
LEM-06 does NOT close by any weighted-energy/Schur/flux method (proven dead, not just
undone). The leap is rescoped to a non-self-adjoint **spectral/Fredholm** bound (resolvent
or Mellin/Wiener–Hopf symbol at `r=−(λ+1)/2`). Circularity firewall: `H(g)(1)=0` is
automatic for L²_W modes — must NOT be asserted from the benign spectrum (NS-C6/SQG trap).

**State as of 2026-11-21 (after independent reviews of the manuscript).** The ledger rows above were last
edited on 2026-07-16 and stay as they are; this paragraph records what has changed around them.

- The manuscript (`paper/stability_spectra_singular_houluo.tex`) now states two theorems: a resolvent bound on
  vertical lines plus point-spectrum uniqueness for the linearization on the constrained domain `Ĝ(½)=0`
  (W1, W2), and small-data support-interior (`X > 1`) nonlinear stability in the edge-analytic class under Gauge (a) with `h ≡ 0` (leaving the distributional `Ω`-equation edge defect `(π/2)(B(τ)-1)δ(X-1)` at `X = 1` whenever `B(τ) ≠ 1`; manuscript Remark 2).
  In this ledger's taxonomy both are at most `proved-on-paper`: the written proofs are in the manuscript's
  appendices, they are not machine-checked, and they have not been independently verified (open items:
  the implicit domain of `M_σ`, the W2 gluing step, the local well-posedness appendix, the nonlinear
  constants). **No per-claim ledger entries exist for them yet.**
- **W3 (semigroup generation + Gearhart–Prüss) is withdrawn** — `docs/desk/w3_gearhart_pruss.md`, Correction.
  The LEM-06 row's "remaining = … W3 generation+GP" is therefore still open, and LEM-06 stays `conjectured`.
- The Lean files added on 2026-11-20/21 (`SlitSchurMargin`, `LinearStrictGapClosure`, `NonlinearComoving`,
  `NonlinearBootstrap`) verify real-number and integer arithmetic and one ODE lemma. Operator bounds enter as
  hypotheses, so none of them moves an analytic claim to `verified`.
- OBS-09 (new) limits what LEM-08 and the firewall sentence above cover: both concern the frozen-exponent
  operator on its constrained domain, not the Gauge-(a) operator.
- THR-01 remains `speculative`. Conjecture 2.4 as stated by Chen–Huang–Li concerns degenerate smooth initial
  data, which the manuscript's perturbation class does not contain.
- **Adversarial review of `0dbf217` (2026-11-21):** `docs/findings_session_2026_11_21_adversarial_review.md`, probe
  `validate/adversarial_review_2026_11_21_probe.py`. Relevant to this ledger: `{H[g](1) = 0}` is not invariant
  (`d/dτ H[g](1) = −H[g](1) − c₁/3`), so the weak edge relation `B = A²` is not kept by any nontrivial solution in the
  manuscript's ansatz at linear order; and the weighted-norm transient of `M_GG` grows with edge depth (peak 12.7 at
  γ = 0.75, `z_min = −14`), which the OBS rows above do not record. No ledger status changes.
