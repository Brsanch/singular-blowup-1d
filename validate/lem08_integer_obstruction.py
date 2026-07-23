#!/usr/bin/env python3
"""
GATE: LEM-08 integer obstruction (discrete core + Stirling budget).

Exercises the shipped module src/lem08_integer_obstruction.py — not a reimplementation.

Checks:
  I1  Exact budget 3π/2 / 2π = 3/4
  I2  Only (k_max,k_min)=(0,0) is L²-admissible among integers in [-4,4]
  I3  integer_obstruction_holds() is True
  I4  Edge-kill with d=0 is impossible; with d≥1 is non-L²
  I5  Measured |Ĝ_p| line decay ≥ 0.92 · (3π/2) at lam=0.7 (feeds the budget)

Exit 0 iff all PASS.
"""
from __future__ import annotations

import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from src.lem08_integer_obstruction import (  # noqa: E402
    L2_INTEGER_BUDGET,
    LAURENT_GROWTH_RATE,
    STIRLING_DECAY_RATE,
    admits_l2_line,
    integer_obstruction_holds,
    laurent_degree,
    l2_integer_budget,
    only_l2_admissible_mode_is_constant,
    point_eigenfunction_impossible,
    stirling_budget_verified,
)

FAILURES: list[str] = []


def report(name: str, ok: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}  {detail}")
    if not ok:
        FAILURES.append(name)


print("I1: exact L² integer budget = 3π/2 ÷ 2π = 3/4")
budget = l2_integer_budget()
ok_budget = (
    abs(budget - 0.75) < 1e-15
    and abs(STIRLING_DECAY_RATE - 1.5 * math.pi) < 1e-15
    and abs(LAURENT_GROWTH_RATE - 2.0 * math.pi) < 1e-15
    and abs(L2_INTEGER_BUDGET - STIRLING_DECAY_RATE / LAURENT_GROWTH_RATE) < 1e-15
)
report("I1 budget identity", ok_budget, f"budget={budget}")

print("I2: only (0,0) is L²-admissible among integer modes")
ok_only = only_l2_admissible_mode_is_constant(range(-4, 5))
report(
    "I2 only constant mode L2",
    ok_only,
    f"admits(0,0)={admits_l2_line(0, 0)} admits(1,0)={admits_l2_line(1, 0)}",
)

print("I3: full discrete obstruction")
ok_full = integer_obstruction_holds()
report("I3 integer_obstruction_holds", ok_full)

print("I4: edge-kill incompatibility cases")
# d=0, try kill edge → impossible
ok_d0 = point_eigenfunction_impossible(0, 0, kills_edge_const=True)
# d=1 modes not L2
ok_d1a = (not admits_l2_line(1, 0)) and point_eigenfunction_impossible(
    1, 0, kills_edge_const=True
)
ok_d1b = (not admits_l2_line(0, -1)) and point_eigenfunction_impossible(
    0, -1, kills_edge_const=True
)
ok_deg = laurent_degree(1, -1) == 2
report(
    "I4 edge-kill cases",
    ok_d0 and ok_d1a and ok_d1b and ok_deg,
    f"d0_imp={ok_d0} d1+={ok_d1a} d1-={ok_d1b} deg(1,-1)={laurent_degree(1, -1)}",
)

print("I5: measured Stirling decay of |Ĝ_p| vs theory 3π/2")
ok_st, rate = stirling_budget_verified(lam=0.7, min_fraction_of_theory=0.92)
report(
    "I5 Stirling line decay",
    ok_st,
    f"measured={rate:.6f} theory={STIRLING_DECAY_RATE:.6f} "
    f"frac={rate / STIRLING_DECAY_RATE:.4f}",
)

print("=" * 72)
if FAILURES:
    print(f"GATE FAIL: {len(FAILURES)} check(s): {FAILURES}")
    sys.exit(1)
print("GATE PASS: LEM-08 integer obstruction (I1–I5)")
print(
    f"KEY: budget={budget} measured_decay={rate:.6f} "
    f"theory_3pi/2={STIRLING_DECAY_RATE:.6f} obstruction={ok_full}"
)
sys.exit(0)
