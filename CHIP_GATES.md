# CHIP_GATES — singular-blowup-1d Lean track

Read before any Lean chip in `lean/`. Ported from the sqg/ns/jacobian discipline
(anti-paraphrase + vacuity gate). The Lean track formalizes the **closed algebraic
backbone** of the Conjecture-2.4 stability proof (Track 1); the analytic coercivity
estimate (the (★) leap, LEM-06) is **Track 2** and is deliberately NOT a `sorry`.

## Build discipline (CLAUDE.md kernel-panic rules — load-bearing)

This M3 Ultra kernel-panics on un-serialized / un-throttled heavy Lean ops. Route
**every** lake op through the machine-global mutex, run **solo**:

- Per file:  `~/.claude/bin/lakelock lake env lean SingularBlowup/Foo.lean`
- Full graph / merge gate:  `~/.claude/bin/lakelock lake build`  (auto-throttles)
- Cache warm:  `~/.claude/bin/lakelock lake exe cache get`  (solo; kill other builds first)
- Never `du` on `lean/.lake/`. Never two heavy ops at once. Pin: mathlib **v4.29.0**.

## Merge gate (all must pass, on a clean worktree)

1. **Full build green:** `~/.claude/bin/lakelock lake build` (elaborates every file
   fresh — catches cross-file duplicate `namespace.decl` that single-file
   `lake env lean` cannot).
2. **Vacuity lint:** `python3 ../noethersolve/scripts/lean_vacuity_lint.py lean/SingularBlowup --max-findings 0`
   — ratcheted baseline **0**, never increases. Flags the SQG `True`-stub class
   (`Prop := True`, `True.intro`, subject-free hypothesis structs, `_`-unused named
   hypotheses). "Zero sorry/axioms" alone is NOT the bar.
3. **Axiom audit:** `#print axioms <decl>` on every new theorem ⇒ only
   `[propext, Classical.choice, Quot.sound]`. Any `sorryAx` = fail.
4. **Aggregator wiring:** every new file imported into `SingularBlowup.lean`
   (un-imported file = vacuous green; `lake build` won't touch it).

## Anti-paraphrase gates (the frontier, not the bloat)

- A new `def Prop`/`structure` whose only discharge is `True`/`rfl`/on the profile
  alone is a renamed sorry. The backbone lemmas (LEM-01..05) are *real theorems*
  about explicit functions — keep that bar.
- "From N inputs" / "under hypothesis H" reformulations of THR-01 (Conjecture 2.4)
  are paraphrases; they do not flip its `speculative` status.
- The actual Track-2 frontier (`docs/ledger/LEM-06`): the (★) bulk-coercivity
  kernel estimate + the H¹_W coupling closure. Formalizing it needs the line
  Hilbert transform / weighted-Sobolev / C₀-semigroup infra mathlib lacks — do NOT
  `sorry`-stub it and call it progress. mathlib-first: grep `.lake/packages/mathlib/`
  before writing any analysis infrastructure.

## Output discipline

Record results as ledger entries (`docs/ledger/`, schema in
`../references/proof_obligation_ledger_schema.md`), not by growing prose docs.
THR-01 caps at `speculative`; LEM-06 stays `conjectured` with OBS gates as evidence;
only per-entry Lean-verified zero-sorry/zero-axiom claims earn `verified`.
