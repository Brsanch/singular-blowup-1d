# NEXT SESSION — complete the rigorous proof + Lean-check it

Goal: turn the numerical stability result for the singular Hou–Luo profile (Conjecture 2.4)
into a rigorous proof, and formalize/verify the checkable parts in Lean. Read first:
`HANDOFF.md` (state-of-play table), `docs/rigor_energy_estimate_setup.md` (the inequalities),
`docs/findings_conj24_evidence_and_frontier.md` (the evidence chain),
`docs/findings_scoop_check_and_rigor_route.md` (route = analytical finite-codimension, NOT a
spectral CAP; template Chen arXiv:2605.15130). The draft writeup of the (done) numerical
contribution is `paper/stability_spectra_singular_houluo.md`.

What is already CLOSED (don't redo): operator form `M_GG=−2∂_z+½Q_op` (verified 1e-11);
energy identity; transport dissipation via exact IBP (`−∫φ'(z)W G²`, inflow boundary
vanishes by Robin, outflow ≤0); edge regularity of the forcing; neutral subspace
`𝒩={ξ@1,µ,ν}`; transient growth (contraction, `µ_an<0`, resolution-independent); and the
**bulk coercivity (★) numerically certified** (bulk energy abscissa `<0` with an explicit
bounded-at-edge anisotropic weight, `validate/energy_estimate_search.py`).

---

## PART A — complete the rigorous proof (the analytic leap)

The estimate to prove (one line, `rigor_energy_estimate_setup.md` §7): there exist an
anisotropic weight `W=e^{−φ}` (`φ'>0`, `φ'→2γ` in the bulk, `φ'→0` at the edge so `W` is
bounded there) and `δ>0` with, for all `(G,h)∈H¹_W`, `G,h ⊥_W 𝒩`, edge condition
`G(−∞)=h(−∞)=0`:
`⟨(G,h), 𝓛(G,h)⟩_{H¹_W} ≤ −δ‖(G,h)‖²`. Recommended order (each numerically de-risked):

- **L1 — bulk coercivity of (★), rigorous.** Prove `⟨G, ½Q_op G⟩_W ≤ ∫φ'W G² − δ‖G‖²_W` for
  `G` supported in the bulk (`z>z_cut`). Numerics: holds, abscissa `−0.02→−0.41`. Tool: an
  **explicit kernel estimate** on `Q_op G = δU/(Y−1)`, `δU=∫_1^Y H(g)` — NOT the pessimistic
  field-of-values (which diverges; `cap_scaffold`). Use the Hilbert-transform identities and
  the averaging structure of `Q_op` (`Q_op G→½H(g)(1)` at the edge, regular). This is the
  core lemma.
- **L2 — edge-layer matched asymptotic — ATTACKED 2026-06-13, RESCOPED (`docs/desk/edge_layer_estimate.md`).**
  6 independent derivation lenses + 18 adversarial verifications + independent numerical
  reproduction (`validate/edge_layer_obstruction.py`, OBS-03). VERDICT: the edge layer does
  **NOT** close by any weighted-energy / Schur / flux method — *proven dead*, not merely undone.
  The nonlocal forcing inherits the base-profile `(Y−1)^{−1/2}` amplitude at the edge
  (`A_full·e^{z_min/2}=0.099` const); the only weight that would dominate re-triggers the
  convective field-of-values divergence (cap_scaffold). The `+40` is a continuum non-normal
  resonance on the g-unbounded `(Y−1)^{−1/2}` branch (κ=1e-3 g-penalty crushes it to negative),
  excluded by admissibility but only numerically. **The leap is rescoped from an energy estimate
  to a non-self-adjoint SPECTRAL/FREDHOLM bound:** "the L²_W realization of `M_GG ⊖ 𝒩` has no
  `Re λ≥0` spectrum", via a verified resolvent bound OR the **Mellin/Wiener–Hopf symbol of `½Q_op`
  at `r=−(λ+1)/2`** (using the continuum identity `H((Y−1)^{−1/2})=0` the discrete operator
  violates). ⚠️ Circularity firewall: `H(g)(1)=0` is automatic for L²_W modes — do NOT assert it
  from the benign spectrum (the NS-C6/SQG circular move). Full route table + next computations in
  the desk doc.
- **L3 — `H¹_W` closure of the coupling.** The `L²` energy does NOT close: `M_Gh=e^{−z/2}∂_z`
  needs `∂_z h` control. Differentiate the system, re-run L1/L2 on `(∂_z G,∂_z h)`, control
  the commutator `[∂_z,F]` and the `e^{−z/2}` factor. This is Chen's "anisotropic estimates +
  IBP exploiting the equation twice" (arXiv:2605.15130) — read its stability section.
- **L4 — nonlinear step.** Linear contraction (L1–L3) + the quadratic self-advection
  `−δU·δΩ_Y` ⇒ nonlinear asymptotic stability modulo `𝒩` (finite-codimension), via a standard
  bootstrap. Numerics: the co-evolving march already relaxes at finite amplitude.

Method discipline (repo): single-agent serial; verify every quantitative constant against the
`validate/` gates as you go; do NOT label anything "proven" until the analytic step is
airtight (this codebase's history: NS Seregin circularity, SQG hypotheses that were `True`).
Each lemma should reduce a `validate/` number to an inequality with an explicit constant.

---

## PART B — Lean-check the proof (honest scope first)

> ✅ **TRACK 1 DONE (2026-06-13).** All five backbone targets are formalized,
> sorry-free, axiom-audited (`[propext, Classical.choice, Quot.sound]` only), and
> pass the full `lake build` + vacuity lint (0 findings). Project: `lean/`, pin
> mathlib **v4.29.0**. Lemmas: `lean/SingularBlowup/Profile.lean`
> (steady residual L1, edge balance L2, c_ω-absorption L3) and
> `lean/SingularBlowup/Energy.lean` (weighted IBP L4, indicial roots L5). Ledger:
> `docs/ledger/` (LEM-01..05 `verified`; LEM-06 (★) `conjectured`; THR-01
> `speculative`). Gates: `CHIP_GATES.md`. Findings:
> `docs/findings_lean_track1_backbone.md`.
> **Next frontier = Track 2 / the analytic leap** = LEM-06: the (★) bulk-coercivity
> kernel estimate + the H¹_W coupling closure (PART A, L1–L4). That is NOT a Lean
> chip yet (mathlib lacks the Hilbert/weighted-Sobolev/semigroup infra) — it is the
> analytic proof, with the `validate/` gates as the numerical guide.

⚠️ **Reality check before writing any Lean.** A full Lean formalization of L1–L4 needs
infrastructure mathlib does **not** currently have in usable form: the line Hilbert transform
as a bounded operator with its identities, weighted Sobolev/energy spaces, and `C₀`-semigroup
/ spectral-gap theory for non-self-adjoint generators. Building that is genuine (in-scope,
per `feedback_mathlib_work_is_in_scope`) but large — do **not** try to Lean the continuous
estimate first, and do **not** create `sorry`-stubbed theorems for the hard parts and call it
progress (`feedback_no_stub_green`). Mathlib-first: grep `.lake/packages/mathlib/` before
writing any analysis infrastructure.

**Two-track Lean strategy (do Track 1 first):**

**Track 1 — the algebraic / finite backbone (genuinely Lean-tractable now).** Self-contained
real-analysis facts that need no Hilbert/Sobolev machinery. Concrete first targets (state and
PROVE these, in order — each is real, not a stub):
1. `Ω̄, Θ̄, Ū` (closed forms) satisfy the steady eq. (2.4) on `Y>1` with `c_l=2, c_ω=−1`
   (pure algebra/calculus on explicit functions; mirrors `verify_singular_steady.py`).
2. `U(1)=−2` and `c_l+2c_ω=0` (the edge/Θ-balance algebra).
3. The `c_ω`-absorption: in `G=(Y−1)^{1/2}g, z=ln(Y−1)`, the transport `−ā g_Y` plus `c_ω g`
   equals `−2∂_z G` (the `+G/−G` cancellation; verified 1e-11 numerically — make it exact).
4. Transport dissipation identity: `∫ e^{−φ} G(−2∂_z G) dz = −∫φ'e^{−φ}G² − [e^{−φ}G²]` (IBP;
   a `MeasureTheory.integral_…`/`intervalIntegral.integral_deriv` exercise).
5. Edge indicial roots: `g~(Y−1)^r` solves the homogeneous edge ODE iff `r=−(λ+1)/2`; weighted
   admissibility `r≥−1/2 ⟺ Re λ≤0`.
These Lean-check the proof's SKELETON (the parts already closed analytically) and are a clean,
finite deliverable.

**Track 2 — the analytic estimate / CAP (research-grade; scope, don't rush).** Two sub-options,
neither near-term:
- (2a) Formalize L1–L4 directly — needs the mathlib Hilbert/weighted-Sobolev/semigroup infra
  above. Long horizon.
- (2b) **Verified-interval-arithmetic CAP of the *discretized* coercivity** (the actual
  inequality `λ_max(−diag(φ')+sym_W(½Q_op_h))<0` on a finite grid), Lean-checked via interval
  arithmetic, plus an analytic truncation/tail bound (the tail bound stays on paper). This is
  how CAP-of-blowup proofs are verified in practice. Check what interval-arithmetic / verified
  linear algebra exists in Lean/mathlib first; this is the more realistic rigorous-Lean route
  but still a project.

**Worked reference architecture for the CAP/spectrum leg** (low-vis harvest deep-read
2026-06-11/12; full notes `../ns-lowvis-harvest/findings/lowvis-taklab-org-RigorousIntegrator-higher-dim-PDEs-20260611.md`
and `...-GC-ubs-CAP-20260611.md`): `taklab-org/RigorousIntegrator_higher_dimensional_PDEs`
(arXiv:2402.00406, J. Sci. Comput. 2024, MATLAB+INTLAB) contains a worked stack for exactly
the "finitely many leading eigenpairs + spectral-gap remainder" architecture: Newton–Kantorovich
validated eigenpairs in weighted ℓ¹ (`Code_Gab/Code Steady State and Eigenvalue Problem
Specific Fm/F_eigenpair.m` etc.) + **Gershgorin-after-similarity** + tail bound for the rest of
the spectrum (`Gershgorin_ring.m`). Gershgorin-after-similarity is a non-normality-robust
enclosure — i.e. the eigenvalue-enclosure route `cap_scaffold` concluded is mandatory after the
numerical-range/energy route diverged. Companion `taklab-org/GC-ubs-CAP` (arXiv:2103.12390) is
the certificate-pattern reference for blow-up CAP (stable-manifold radii polynomials +
validated-Lyapunov-domain absorption + blow-up-time integral). What does NOT transfer = the
actual new work: their bases are cosine-on-rectangles / Taylor-at-equilibria (we need the
half-line `z`-grid with algebraic edge behavior in the `e^{−2γz}` weight), and their tail bound
is sectorial-semigroup (`M_GG=−2∂_z+½Q_op` is transport-dominated, non-sectorial). ⚠️ Port WITH
tests only: our deep-read found 9 rigor defects in the shipped code (sign error in a
radii-polynomial Z2 tail; trapping radius computed with float `roots` on `mid()`-stripped
intervals; unchecked `abs()` that would silently accept an unstable rate as decay) — re-derive
every hand-assembled bound.

**Project setup steps (next session):**
- `mkdir lean && cd lean && lake new SingularBlowup` (or copy the lakefile/toolchain pin from
  `../sqg-lean-proofs` to match a known-good mathlib version — check its `lean-toolchain` and
  `lakefile`).
- ⚠️ **Build discipline (CLAUDE.md kernel-panic rules):** never run unthrottled `lake build`
  or `lake exe cache get`. Per-file: `LEAN_NUM_THREADS=1 lake env lean FILE.lean`. Full graph /
  merge gate: `taskpolicy -b nice -n 19 env LEAN_NUM_THREADS=1 lake build`.
- Put Track-1 lemmas in `SingularBlowup/Profile.lean` (steady state, U(1), c_ω-absorption) and
  `SingularBlowup/Energy.lean` (IBP identity, indicial roots). Keep files small (1–3k LOC,
  textbook granularity), not 100-LOC chips.
- Adopt the proof-obligation ledger schema from day one
  (`../references/proof_obligation_ledger_schema.md`, ported 2026-06-12 from the harvest):
  Track-1 lemmas = LEM entries with `machine_verification: Lean`; the (★) estimate and L1–L4
  stay `conjectured` (with the `validate/` gates as OBS evidence) until proved; the headline
  stability claim stays `speculative` + confidence until external review.

**Honest Lean verdict:** Track 1 is a real, bounded, achievable deliverable that formalizes the
closed skeleton. Track 2 (the actual coercivity/the full theorem) is research-grade and
gated on mathlib analysis infrastructure or a Lean interval-arithmetic CAP — set expectations
accordingly; do not present a `sorry`-d full theorem as a result.

---

## Reading order for the next session
1. `HANDOFF.md` (table) → `docs/rigor_energy_estimate_setup.md` (the inequalities, §4 the bulk
   result, §7 the one-line target).
2. `docs/findings_cap_scaffold.md` + `findings_transient_growth_edge_conditioned.md` (why the
   energy route works in the bulk; the log-norm = energy abscissa).
3. Chen arXiv:2605.15130 stability section (the finite-codimension + anisotropic + IBP-twice
   template) — fetch the full PDF for the method.
4. `validate/energy_estimate_search.py`, `coevolve_*.py` (the numerical de-risking + constants).

Numerical contribution (done) → write up from `paper/stability_spectra_singular_houluo.md`.
