"""
LEM-08 integer obstruction — pure, unit-testable core.

The point-spectrum rigidity argument (docs/desk/edge_mellin_symbol.md §5.3–5.5,
validate/edge_pointspectrum_rigid.py) reduces to exact integers once the Mellin
line-decay of the particular solution Ĝ_p is known:

  |Ĝ_p(c + iτ)| ∼ exp(−(3π/2)|τ|)   (λ-independent leading exponential)

A genuine 1-periodic entire factor P is a Laurent polynomial in w = exp(2π i s)
of degree d = k_max − k_min. On the vertical line it grows as
  exp(2π k_max |τ|) or exp(2π |k_min| |τ|)
on at least one side when d ≥ 1. Line-L² of Ĝ_p · P therefore requires
  2π k_max < 3π/2  and  2π |k_min| < 3π/2
i.e.  k_max < 3/4  and  |k_min| < 3/4.
For integers this forces k_max = k_min = 0, hence d = 0 (no zeros).

But killing the s = 0 edge-constant pole (residue ∝ Γ(λ) ≠ 0 for Re λ ≥ 0,
λ ∉ {0,−1,…}) needs at least one zero per period strip ⇒ d ≥ 1.

Incompatibility ⇒ no L²-line / L²_W eigenfunction for Re λ ≥ 0 except the
integer-λ = 𝒩 termination cases (quotiented out).

This module encodes only the **integer / budget half** (no Gamma construction).
The decay constant is verified numerically in the companion gate against Ĝ_p.
"""
from __future__ import annotations

import math
from typing import Iterable

# Half-plane decay budget of |Ĝ_p| on a vertical Mellin line (leading Stirling).
# Units: |Ĝ_p(c+iτ)| ≲ C exp(−STIRLING_DECAY_RATE |τ|).
STIRLING_DECAY_RATE: float = 1.5 * math.pi  # 3π/2

# Growth per integer Laurent mode on the line: |w^k| = exp(2π k τ) in magnitude.
LAURENT_GROWTH_RATE: float = 2.0 * math.pi  # 2π

# Integer bound: k_max < STIRLING_DECAY_RATE / LAURENT_GROWTH_RATE = 3/4.
L2_INTEGER_BUDGET: float = STIRLING_DECAY_RATE / LAURENT_GROWTH_RATE  # 0.75


def l2_integer_budget() -> float:
    """Return the exact L²-line budget k_max < 3/4 (and |k_min| < 3/4)."""
    return L2_INTEGER_BUDGET


def admits_l2_line(k_max: int, k_min: int) -> bool:
    """
    Whether a Laurent factor with modes k_min..k_max keeps Ĝ_p·P square-integrable
    on a vertical Mellin line, given the 3π/2 decay of Ĝ_p.
    """
    if k_max < k_min:
        raise ValueError("k_max must be ≥ k_min")
    return (k_max < L2_INTEGER_BUDGET) and (abs(k_min) < L2_INTEGER_BUDGET)


def laurent_degree(k_max: int, k_min: int) -> int:
    """Degree d = k_max − k_min of a Laurent polynomial supported on [k_min, k_max]."""
    if k_max < k_min:
        raise ValueError("k_max must be ≥ k_min")
    return k_max - k_min


def kills_edge_const_requires_degree_at_least_one(d: int) -> bool:
    """
    A nonzero 1-periodic entire factor that vanishes at one coset per strip
    (needed to cancel the s=0 edge-constant pole) requires d ≥ 1.
    """
    return d >= 1


def point_eigenfunction_impossible(
    k_max: int,
    k_min: int,
    *,
    kills_edge_const: bool,
) -> bool:
    """
    True iff the (k_max, k_min, edge-kill) triple is incompatible with an
    L²-line / L²_W eigenfunction under the LEM-08 budget.

    Cases:
      - admits L² and does not kill edge → not an eigenfunction (edge mass infinite)
      - does not admit L² → not an eigenfunction (line divergent)
      - admits L² and kills edge → impossible for integers (d=0 vs d≥1)
    """
    if k_max < k_min:
        raise ValueError("k_max must be ≥ k_min")
    d = laurent_degree(k_max, k_min)
    l2_ok = admits_l2_line(k_max, k_min)
    if kills_edge_const:
        # Need d ≥ 1 to place a zero; L² forbids any nonzero mode.
        return (not l2_ok) or (d == 0)
    # Without killing the edge const, either line fails or edge mass is infinite.
    return True


def enumerate_integer_modes(
    k_range: Iterable[int] = range(-2, 3),
) -> list[tuple[int, int, bool, bool, bool]]:
    """
    Exhaustive scan of small integer (k_max, k_min) pairs with k_max ≥ k_min.

    Returns rows (k_max, k_min, l2_ok, d_ge_1, impossible_with_edge_kill).
    """
    ks = list(k_range)
    rows: list[tuple[int, int, bool, bool, bool]] = []
    for k_max in ks:
        for k_min in ks:
            if k_max < k_min:
                continue
            l2_ok = admits_l2_line(k_max, k_min)
            d = laurent_degree(k_max, k_min)
            d_ge_1 = kills_edge_const_requires_degree_at_least_one(d)
            imp = point_eigenfunction_impossible(
                k_max, k_min, kills_edge_const=True
            )
            rows.append((k_max, k_min, l2_ok, d_ge_1, imp))
    return rows


def only_l2_admissible_mode_is_constant(
    k_range: Iterable[int] = range(-3, 4),
) -> bool:
    """
    Among integer modes in k_range, the only L²-admissible Laurent support
    is (k_max, k_min) = (0, 0).
    """
    admissible = [
        (k_max, k_min)
        for k_max, k_min, l2_ok, _, _ in enumerate_integer_modes(k_range)
        if l2_ok
    ]
    return admissible == [(0, 0)]


def integer_obstruction_holds() -> bool:
    """
    Full discrete LEM-08 obstruction on a neighborhood of the origin:
      - only (0,0) is L²-admissible
      - killing the edge const with (0,0) fails (d=0)
      - every d≥1 pair is not L²-admissible
    ⇒ no compatible eigenfunction mode.
    """
    if not only_l2_admissible_mode_is_constant(range(-4, 5)):
        return False
    # Constant mode cannot kill the edge pole.
    if not point_eigenfunction_impossible(0, 0, kills_edge_const=True):
        return False
    # Every nonzero-degree pair in a window is either non-L² or still impossible.
    for k_max, k_min, l2_ok, d_ge_1, imp in enumerate_integer_modes(range(-3, 4)):
        d = laurent_degree(k_max, k_min)
        if d == 0:
            continue
        if l2_ok:
            return False  # would break the 3/4 budget for integers
        if d_ge_1 and not imp:
            return False
    return True


# ---------------------------------------------------------------------------
# Numerical Stirling budget for Ĝ_p (feeds L2_INTEGER_BUDGET; not pure-integer)
# ---------------------------------------------------------------------------

def ghat_p_abs(s: complex, lam: float) -> float:
    """|Γ(s+λ) Γ(s+2) Γ(s/2) Γ(1−s/2)| via cmath-friendly magnitude product."""
    # Local import keeps the pure-integer path free of scipy for unit tests.
    from scipy.special import gamma as G

    val = G(s + lam) * G(s + 2) * G(s / 2) * G(1 - s / 2)
    return float(abs(val))


def measured_decay_rate(
    lam: float = 0.7,
    c: float = 0.5,
    tau_lo: float = 8.0,
    tau_hi: float = 16.0,
) -> float:
    """
    Estimate the exponential decay rate of |Ĝ_p(c+iτ)| between |τ|=tau_lo and tau_hi:
      rate ≈ −(log|G_hi| − log|G_lo|) / (tau_hi − tau_lo).
    Should approach 3π/2 ≈ 4.712 from below for moderate τ.
    """
    g_lo = ghat_p_abs(complex(c, tau_lo), lam)
    g_hi = ghat_p_abs(complex(c, tau_hi), lam)
    if g_lo <= 0.0 or g_hi <= 0.0:
        raise RuntimeError("Ĝ_p vanished unexpectedly on the sample line")
    return -(math.log(g_hi) - math.log(g_lo)) / (tau_hi - tau_lo)


def stirling_budget_verified(
    *,
    lam: float = 0.7,
    min_fraction_of_theory: float = 0.92,
) -> tuple[bool, float]:
    """
    Check that the measured line decay is at least min_fraction_of_theory * (3π/2).
    Returns (ok, measured_rate).
    """
    rate = measured_decay_rate(lam=lam)
    ok = rate >= min_fraction_of_theory * STIRLING_DECAY_RATE
    return ok, rate
