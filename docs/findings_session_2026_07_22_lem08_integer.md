# Session 2026-07-22 — LEM-08 integer obstruction (machine + Lean)

## Attack choice (5 lines)

1. **Primary target:** LEM-08 point-spectrum integer obstruction (ledger `proved-on-paper`, no Lean yet).
2. **Why:** HANDOFF ordered list puts linear-gap W1–W3 as execution, then **Lean LEM-08** as secondary; W1 is write-up-heavy; LEM-08 has an exact discrete core that is firewall-clean and chip-sized.
3. **First tool:** `validate/edge_pointspectrum_rigid.py` baseline (PASS) + new pure module for the growth budget + degree incompatibility.
4. **Deliverable:** pure functions in `src/lem08_integer_obstruction.py`, gate `validate/lem08_integer_obstruction.py`, Lean `SingularBlowup/PointSpectrum.lean` (integer half only — no Gamma formalization).
5. **Not claiming:** essential-spectrum gap, Conjecture 2.4, semigroup/Gearhart–Prüss, or nonlinear stability.

## Scope honesty

This closes the **discrete incompatibility** (line-`L²` ⇒ `d=0`; edge-const kill ⇒ `d≥1`) and machine-checks the **Stirling line-decay budget** `3π/2` that forces `k_max < 3/4`. It does **not** Lean-prove the Gamma asymptotics or the full eigenfunction construction (those remain paper + mpmath gates).

## Results (machine-checkable)

### Python gate `validate/lem08_integer_obstruction.py` (2× PASS)

| Check | Claim | Measured |
|-------|--------|----------|
| I1 | budget `3π/2 ÷ 2π = 3/4` | `budget=0.75` |
| I2 | only `(0,0)` L²-admissible | `admits(0,0)=True`, `admits(1,0)=False` |
| I3 | `integer_obstruction_holds()` | `True` |
| I4 | edge-kill incompat | d=0 impossible; d≥1 non-L² |
| I5 | Stirling decay ≥ 0.92·(3π/2) | **measured=4.480288**, theory=4.712389, frac=**0.9507** |

Shipped entry point: `src/lem08_integer_obstruction.py` (pure budget + obstruction; gate imports it — no reimplementation).

### Lean `SingularBlowup/PointSpectrum.lean` (Track-1 extension)

- Theorems: `int_lt_three_quarters_eq_zero`, `only_constant_mode_l2_admissible`, `integer_obstruction`, `no_l2_and_edge_kill`, `constant_mode_cannot_kill_edge`.
- Aggregator: imported from `SingularBlowup.lean`.
- `lake env lean SingularBlowup/PointSpectrum.lean` → **exit 0** (throttled, LEAN_NUM_THREADS=1).
- Vacuity lint on `lean/SingularBlowup`: **0 findings**.

### Not claimed

- Essential-spectrum strict gap, W1 Plemelj write-up, Lemma E′, Gearhart–Prüss, nonlinear/co-moving.
- Full LEM-08 is still `proved-on-paper` for the continuous half; Lean now covers the **integer skeleton** (`lean_ref` updated).

## Next (ordered)

1. W1 Plemelj/Stirling write-up (desk §10)
2. Lemma E′ edge bootstrap
3. W3 semigroup + Gearhart–Prüss
