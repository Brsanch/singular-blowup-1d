"""
Quantitative dressed two-weight Schur envelopes for the slit-plane resolvent
(docs/desk/convective_weight_spectrum.md §10, W1 constants).

Pure functions: (sigma, delta, tau-grid) -> envelope scalars / tables.
delta = Re(lambda) + 2*sigma > 0 is the distance of lambda above the transport line.

Kernel matches validate/edge_slitplane_representation.py:
  dressed(tau, tau') = |X(s)/X(s')| * ker(tau,tau') / |lambda + 2 s'|
  ker = 4 pi |w'|/|w'-w|  (plain on tau'>0 wing; * |w|/|w'| on tau'<0 subtraction wing)
  w = exp(4 pi i s), s = sigma + i tau

No discrete Hilbert; Gamma magnitudes only (Stirling-safe via cmath/mpmath).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np

try:
    import mpmath as mp
except ImportError:  # pragma: no cover
    mp = None


def _gamma_abs_mp(z: complex, dps: int = 25) -> float:
    """|Gamma(z)| via mpmath (stable for moderate |Im z|)."""
    assert mp is not None
    with mp.workdps(dps):
        return float(abs(mp.gamma(mp.mpc(z.real, z.imag))))


def X_abs(s: complex, lam: complex, dps: int = 25) -> float:
    """|X(s)| with X(s) = Gamma(2s+lam) Gamma(2s+2) Gamma(s) Gamma(1-s)."""
    return (
        _gamma_abs_mp(2 * s + lam, dps)
        * _gamma_abs_mp(2 * s + 2, dps)
        * _gamma_abs_mp(s, dps)
        * _gamma_abs_mp(1 - s, dps)
    )


def wray_abs_ratio(tau: float, taup: float, sigma: float) -> float:
    """
    Split-kernel factor ker(tau, taup) = 4 pi |w'| / |w' - w|
    times subtraction |w|/|w'| when taup < 0.
    w = exp(4 pi i (sigma + i tau)) = exp(4 pi i sigma) exp(-4 pi tau).
    So |w| = exp(-4 pi tau), |w'| = exp(-4 pi taup).
    """
    # |w' - w|^2 with common phase factored out:
    # w' - w = e^{4 pi i sigma} ( e^{-4 pi taup} - e^{-4 pi tau} )
    # |w' - w| = |e^{-4 pi taup} - e^{-4 pi tau}|
    a = math.exp(-4 * math.pi * taup)
    b = math.exp(-4 * math.pi * tau)
    # numerically stable when both large/small
    if taup >= tau:
        # a <= b; factor a: |1 - exp(-4pi (tau-taup))|
        diff = a * abs(1.0 - math.exp(-4 * math.pi * (tau - taup)))
    else:
        diff = b * abs(1.0 - math.exp(-4 * math.pi * (taup - tau)))
    if diff == 0.0:
        return float("inf")
    ker = 4 * math.pi * a / diff
    if taup < 0:
        # subtraction wing: ker *= |w|/|w'| = exp(-4 pi (tau - taup))
        ker *= math.exp(-4 * math.pi * (tau - taup))
    return ker


def dressed_kernel(
    tau: float,
    taup: float,
    sigma: float,
    lam: complex,
    dps: int = 25,
    x_cache: Optional[dict] = None,
) -> float:
    """Dressed two-weight Schur kernel magnitude at (tau, taup)."""
    s = complex(sigma, tau)
    sp = complex(sigma, taup)

    def xabs(z: complex) -> float:
        key = (round(z.real, 12), round(z.imag, 12), round(lam.real, 12), round(lam.imag, 12))
        if x_cache is not None and key in x_cache:
            return x_cache[key]
        val = X_abs(z, lam, dps=dps)
        if x_cache is not None:
            x_cache[key] = val
        return val

    ratio = xabs(s) / xabs(sp)
    ker = wray_abs_ratio(tau, taup, sigma)
    denom = abs(lam + 2 * sp)
    if denom == 0.0:
        return float("inf")
    return ratio * ker / denom


def default_tau_grid(tmax: float = 8.0, n: int = 33) -> np.ndarray:
    """Symmetric log-ish grid in [-tmax, tmax] including 0."""
    # denser near 0 for diagonal band; extend to tmax
    pos = np.unique(
        np.concatenate(
            [
                np.linspace(0.0, 1.0, max(5, n // 4)),
                np.linspace(1.0, 3.0, max(5, n // 4)),
                np.linspace(3.0, tmax, max(5, n // 3)),
            ]
        )
    )
    return np.unique(np.concatenate([-pos[::-1], pos]))


@dataclass
class EnvelopeRow:
    sigma: float
    delta: float  # Re(lambda) + 2*sigma
    lam_real: float
    # sup_{|Δτ| >= band} dressed * exp(+pi |Δτ|)
    sup_off_pi: float
    # measured off-diagonal decay rate α from log-dressed vs |Δτ| (large sep)
    alpha_eff: float
    # max dressed on |Δτ| <= band
    max_diag: float
    # max |X(σ+iτ)/X(σ+iτ')| on |Δτ|<=band (weight-ratio part)
    max_weight_ratio_diag: float
    # bare defect density growth rate γ (k ~ e^{+γ|τ|}); Stirling predicts ≈ 3π
    defect_growth: float
    # k(τ=8) scale check
    defect_k8: float
    # termwise ladder sum proxy: sum_n sup_τ |Π_n| * delta  (~1.15 at bulk)
    termwise_C: float
    n_pairs_off: int
    n_pairs_diag: int


def measure_alpha_eff(
    sigma: float,
    lam: complex,
    sep_list: Sequence[float] = (2.0, 4.0, 6.0, 8.0, 10.0),
    dps: int = 25,
) -> float:
    """
    Effective off-diagonal decay on the former-corner wing (both τ,τ' < 0,
    |τ'| > |τ| — the slowest desk-§10 quadrant).

    Fix τ' = -12 and vary τ = τ' + sep. Theoretical Schur margin α ≥ π.
    """
    x_cache: dict = {}
    tp = -12.0
    vals = []
    for sep in sep_list:
        t = tp + float(sep)
        v = dressed_kernel(t, tp, sigma, lam, dps=dps, x_cache=x_cache)
        vals.append((float(sep), v))
    slopes = []
    for i in range(len(vals) - 1):
        d0, v0 = vals[i]
        d1, v1 = vals[i + 1]
        if v0 <= 0 or v1 <= 0 or not math.isfinite(v0) or not math.isfinite(v1):
            continue
        alpha = -(math.log(v1) - math.log(v0)) / (d1 - d0)
        slopes.append(alpha)
    if not slopes:
        return float("nan")
    return float(sum(slopes) / len(slopes))

def termwise_ladder_C(
    sigma: float,
    lam: complex,
    n_max: int = 40,
    tau_grid: Optional[np.ndarray] = None,
    dps: int = 25,
) -> float:
    """
    C_term ≈ (sum_{n=0}^{n_max} sup_τ |Π_n(σ+iτ)|) * dist
    with dist = Re(lam) + 2*sigma, Π_n = X(s)/[X(s+n/2)(lam+2s+n)].
    Desk §8 reports ≈ 1.15 at one (σ,λ); here as function of delta.
    """
    if tau_grid is None:
        tau_grid = default_tau_grid(tmax=6.0, n=25)
    dist = lam.real + 2 * sigma
    if dist <= 0:
        return float("inf")
    x_cache: dict = {}
    total = 0.0
    for n in range(n_max + 1):
        m_sup = 0.0
        for tau in tau_grid:
            s = complex(sigma, float(tau))
            sn = s + 0.5 * n
            key_s = (round(s.real, 12), round(s.imag, 12), round(lam.real, 12), round(lam.imag, 12))
            key_sn = (round(sn.real, 12), round(sn.imag, 12), round(lam.real, 12), round(lam.imag, 12))
            if key_s not in x_cache:
                x_cache[key_s] = X_abs(s, lam, dps=dps)
            if key_sn not in x_cache:
                x_cache[key_sn] = X_abs(sn, lam, dps=dps)
            denom = abs(lam + 2 * s + n)
            if denom == 0:
                pi_n = float("inf")
            else:
                pi_n = x_cache[key_s] / (x_cache[key_sn] * denom)
            if pi_n > m_sup:
                m_sup = pi_n
        total += m_sup
    return total * dist


def defect_density_abs(taup: float, sigma: float, lam: complex, dps: int = 25) -> float:
    """
    Line density for the rank-one defect / Cauchy data (no Fhat):
      k(τ') = 1 / (|X(s')| |λ + 2 s'|) ,  s' = σ + i τ'.
    """
    sp = complex(sigma, taup)
    denom = X_abs(sp, lam, dps=dps) * abs(lam + 2 * sp)
    if denom == 0.0:
        return float("inf")
    return 1.0 / denom


def defect_growth_rate(
    sigma: float,
    lam: complex,
    tau_samples: Sequence[float] = (4.0, 6.0, 8.0, 10.0, 12.0),
    dps: int = 25,
) -> Tuple[float, float]:
    """
    Log-growth rate of bare defect density k(τ)=1/(|X||λ+2s|) on the +τ wing:
        k(τ) ~ C e^{+γ |τ|}  ⇒  γ = + d log k / d|τ|   (report γ; γ>0 means growth).

    Desk §10 prose claimed ~ e^{-3π|τ|/2} *decay*. Stirling on X gives
    |X| ~ e^{-3π|τ|} ⇒ k ~ e^{+3π|τ|} growth (γ ≈ 3π ≈ 9.42). Returns
    (gamma_growth, k_at_tau=8) as a scale check.
    """
    pos = []
    for t in tau_samples:
        t = float(t)
        pos.append((t, defect_density_abs(t, sigma, lam, dps=dps)))
    slopes = []
    for i in range(len(pos) - 1):
        t0, k0 = pos[i]
        t1, k1 = pos[i + 1]
        if k0 <= 0 or k1 <= 0 or not math.isfinite(k0) or not math.isfinite(k1):
            continue
        # growth rate γ with k ~ e^{γ t}
        slopes.append((math.log(k1) - math.log(k0)) / (t1 - t0))
    gamma = float(sum(slopes) / len(slopes)) if slopes else float("nan")
    k8 = defect_density_abs(8.0, sigma, lam, dps=dps)
    return gamma, k8

def compute_envelope_row(
    sigma: float,
    delta: float,
    band: float = 1.0,
    tau_grid: Optional[np.ndarray] = None,
    dps: int = 25,
    n_max_termwise: int = 30,
) -> EnvelopeRow:
    """
    Full envelope row at (sigma, delta) with real lambda = delta - 2*sigma.
    """
    if delta <= 0:
        raise ValueError("delta = Re(lambda)+2*sigma must be > 0")
    if not (0 < sigma < 0.5):
        raise ValueError("sigma must be in (0, 1/2)")
    lam = complex(delta - 2 * sigma, 0.0)
    if tau_grid is None:
        tau_grid = default_tau_grid(tmax=8.0, n=33)
    x_cache: dict = {}
    sup_off = 0.0
    max_diag = 0.0
    max_wr = 0.0
    n_off = 0
    n_diag = 0
    pi = math.pi
    for t in tau_grid:
        for tp in tau_grid:
            if t == tp:
                # diagonal singular in ker; skip exact coincidence
                continue
            dabs = abs(float(t) - float(tp))
            val = dressed_kernel(float(t), float(tp), sigma, lam, dps=dps, x_cache=x_cache)
            if not math.isfinite(val):
                continue
            if dabs <= band:
                n_diag += 1
                if val > max_diag:
                    max_diag = val
                # weight ratio only
                s = complex(sigma, float(t))
                sp = complex(sigma, float(tp))
                wr = X_abs(s, lam, dps=dps) / X_abs(sp, lam, dps=dps)
                if wr > max_wr:
                    max_wr = wr
            else:
                n_off += 1
                env = val * math.exp(pi * dabs)
                if env > sup_off:
                    sup_off = env
    alpha = measure_alpha_eff(sigma, lam, dps=dps)
    def_g, def_k8 = defect_growth_rate(sigma, lam, dps=dps)
    c_term = termwise_ladder_C(sigma, lam, n_max=n_max_termwise, tau_grid=tau_grid, dps=dps)
    return EnvelopeRow(
        sigma=sigma,
        delta=delta,
        lam_real=float(lam.real),
        sup_off_pi=sup_off,
        alpha_eff=alpha,
        max_diag=max_diag,
        max_weight_ratio_diag=max_wr,
        defect_growth=def_g,
        defect_k8=def_k8,
        termwise_C=c_term,
        n_pairs_off=n_off,
        n_pairs_diag=n_diag,
    )


def envelope_table(
    sigmas: Iterable[float] = (0.25, 0.35, 0.4),
    deltas: Iterable[float] = (0.1, 0.3, 0.5, 1.0, 1.5),
    **kwargs,
) -> List[EnvelopeRow]:
    rows = []
    for sigma in sigmas:
        for delta in deltas:
            rows.append(compute_envelope_row(float(sigma), float(delta), **kwargs))
    return rows


def key_numbers(rows: Sequence[EnvelopeRow]) -> dict:
    """
    Compact key numbers for gates / dual-run comparison.
    Key claim: at sigma=0.4, delta=1.1 (lam=0.3, the desk anchor), report
    sup_off_pi, alpha_eff, max_diag, termwise_C, and the delta-slope of
    sup_off_pi near the line.
    """
    by = {(round(r.sigma, 6), round(r.delta, 6)): r for r in rows}
    # pick anchor if present else first
    anchor = by.get((0.4, 1.1)) or by.get((0.4, 1.0)) or rows[0]
    # delta dependence at fixed sigma=0.4: C_off(delta) ≈ a + b/delta?
    s04 = [r for r in rows if abs(r.sigma - 0.4) < 1e-9]
    s04_sorted = sorted(s04, key=lambda r: r.delta)
    # fit sup_off_pi * delta ~ const for small delta
    if len(s04_sorted) >= 2:
        r_small = s04_sorted[0]
        r_mid = s04_sorted[len(s04_sorted) // 2]
        product_small = r_small.sup_off_pi * r_small.delta
        product_mid = r_mid.sup_off_pi * r_mid.delta
    else:
        product_small = anchor.sup_off_pi * anchor.delta
        product_mid = product_small
    return {
        "anchor_sigma": anchor.sigma,
        "anchor_delta": anchor.delta,
        "anchor_sup_off_pi": anchor.sup_off_pi,
        "anchor_alpha_eff": anchor.alpha_eff,
        "anchor_max_diag": anchor.max_diag,
        "anchor_termwise_C": anchor.termwise_C,
        "anchor_defect_growth": anchor.defect_growth,
        "anchor_defect_k8": anchor.defect_k8,
        "sup_off_pi_times_delta_small": product_small,
        "sup_off_pi_times_delta_mid": product_mid,
        "n_rows": len(rows),
    }
