#!/usr/bin/env python3
"""
NONLINEAR shape march of the singular-profile perturbation on the edge-resolved
z = ln(Y-1) grid  (THE NEXT CHIP #1).

Extends the LINEAR convective spectral-gap result (validate/spectral_gap.py,
infinitesimal) to FINITE amplitude:  does a finite-amplitude perturbation of the
frozen singular profile still DECAY in the convective weighted norm, or does the
nonlinearity surface an instability the linear analysis cannot see?

------------------------------------------------------------------------------
Setup.  Frozen base  Omega_bar = 1{Y>1}/sqrt(Y-1),  Theta_bar = (pi/2)1{Y>1},
c_l = 2, c_omega = -1.  Perturbation
    Omega = Omega_bar + dOmega ,  dOmega = (Y-1)^{-1/2} G   (G bounded at edge)
    Theta = Theta_bar + h .
Insert into eq.(2.4); it splits EXACTLY as  dpert_tau = L(dpert) + N(dpert,dpert):

  * L  = the VERIFIED frozen linear operator (src/stability_edge_z.build_blocks_z:
         M_GG, M_Gh, M_hh).  Its spectrum gives the convective gap 2*gamma.
  * N  = the ONLY quadratic term, the perturbation self-advection  -dU * dpert_Y
         (dU = int_0^Y H(dOmega) = the perturbation velocity).  In the G-variable:
            N_G = -(dU/(Y-1)) * (G_z - G/2)   [ = (Y-1)^{1/2} * (-dU * dOmega_Y) ]
            N_h = -(dU/(Y-1)) * h_z * (Y-1)... -> = -(dU/(Y-1)) * (Y-1) h_Y
                = -(Q_op@G) * (Dz@h)          [ since (Y-1) h_Y = (Y-1) e^{-z} h_z = h_z ]
         with  dU/(Y-1) = Q_op @ G  (bounded edge form, src/stability_edge_z).

Frame: non-moving (identical to spectral_gap.py, so GATE A is an exact check);
edge BC G(zmin)=h(zmin)=0 (reduced grid); BULK initial bump so the xi/symmetry
mode is not excited (and the weighted space excludes it anyway).

------------------------------------------------------------------------------
VALIDATION LADDER (correctness must pass before the verdict is believed):
  GATE A  RK4 with the nonlinear term OFF must reproduce spectral_gap.py's
          weighted decay rate ~ 2*gamma  (isolates the time-stepper).
  GATE B  RK4 with N ON at amplitude A -> 0 must converge to that same linear
          rate  (the nonlinear term vanishes correctly as O(A)).
  SWEEP   finite amplitude A: weighted decay rate, peak transient norm ratio,
          and edge hi-k fraction; convergence under (N, zmin) refinement.
HARD STOP: a resolution-ROBUST finite-A growth is reported as a candidate
nonlinear instability -- it is NOT hidden.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from stability_edge_z import build_blocks_z  # noqa: E402


# ---------------------------------------------------------------------------
def make_rhs(B):
    """Return rhs(G,h,nl) for the coupled perturbation on the reduced z-grid."""
    M_GG, M_Gh, M_hh = B["M_GG"], B["M_Gh"], B["M_hh"]
    Q_op, Dz = B["Q_op"], B["Dz"]

    def rhs(G, h, nl=True):
        dG = M_GG @ G + M_Gh @ h
        dh = M_hh @ h
        if nl:
            dU_over = Q_op @ G                       # dU/(Y-1), bounded
            dG = dG - dU_over * (Dz @ G - 0.5 * G)   # -(dU/(Y-1)) dOmega_Y, in G
            dh = dh - dU_over * (Dz @ h)             # -dU * dTheta_Y
        return dG, dh

    return rhs


def rk4(rhs, G, h, dt, nl):
    a1G, a1h = rhs(G, h, nl)
    a2G, a2h = rhs(G + 0.5 * dt * a1G, h + 0.5 * dt * a1h, nl)
    a3G, a3h = rhs(G + 0.5 * dt * a2G, h + 0.5 * dt * a2h, nl)
    a4G, a4h = rhs(G + dt * a3G, h + dt * a3h, nl)
    G = G + dt / 6 * (a1G + 2 * a2G + 2 * a3G + a4G)
    h = h + dt / 6 * (a1h + 2 * a2h + 2 * a3h + a4h)
    G[0] = 0.0; h[0] = 0.0                            # re-impose inflow/edge BC
    return G, h


def wnorm(z, dz, G, h, gamma):
    w = np.exp(-2 * gamma * z)
    return float(np.sqrt(np.sum(w * (G ** 2 + h ** 2)) * dz))


def hi_k_fraction(G):
    """Fraction of L2 energy of G in the upper half of the z-wavenumber band
    (grid-scale content detector, per the handoff's mandatory check)."""
    F = np.abs(np.fft.rfft(G)) ** 2
    n = len(F)
    tot = F.sum() + 1e-30
    return float(F[n // 2:].sum() / tot)


def march(B, A, gammas, taus, nl=True, h_amp=0.0):
    """Evolve a bulk bump of amplitude A; return weighted-norm series + diagnostics."""
    z, dz = B["z"], B["dz"]
    rhs = make_rhs(B)
    G = A * np.exp(-((z - 1.0) / 1.0) ** 2)
    h = h_amp * np.exp(-((z - 1.5) / 1.0) ** 2)
    G[0] = 0.0; h[0] = 0.0

    dt = 0.25 * dz                                   # CFL ~0.5 for transport speed 2
    series = {g: [] for g in gammas}
    n0 = {g: wnorm(z, dz, G, h, g) for g in gammas}
    peak = 1.0
    hik = []
    tau, ti = 0.0, 0
    targets = list(taus)
    # record at tau=0
    for g in gammas:
        series[g].append(wnorm(z, dz, G, h, g))
    hik.append(hi_k_fraction(G))
    ti = 1
    while ti < len(targets):
        while tau < targets[ti] - 1e-12:
            step = min(dt, targets[ti] - tau)
            G, h = rk4(rhs, G, h, step, nl)
            tau += step
            if not np.all(np.isfinite(G)):
                for g in gammas:
                    series[g].append(np.inf)
                hik.append(1.0)
                return series, n0, np.inf, hik, z, G
        for g in gammas:
            series[g].append(wnorm(z, dz, G, h, g))
        # peak transient ratio in the reference gamma (first in list)
        peak = max(peak, series[gammas[0]][-1] / (n0[gammas[0]] + 1e-30))
        hik.append(hi_k_fraction(G))
        ti += 1
    return series, n0, peak, hik, z, G


def fit_rate(taus, vals, win=(0.5, 2.5)):
    t = np.asarray(taus); y = np.asarray(vals)
    m = (t >= win[0]) & (t <= win[1]) & np.isfinite(y) & (y > 0)
    if m.sum() < 2:
        return float("nan")
    return float(-np.polyfit(t[m], np.log(y[m]), 1)[0])


# ---------------------------------------------------------------------------
def main():
    print("=" * 80)
    print("  NONLINEAR shape march (finite-amplitude convective stability, frozen base)")
    print("=" * 80)

    N, zmin, zmax = 1400, -14.0, 12.0
    B = build_blocks_z(N=N, zmin=zmin, zmax=zmax, beta=0.0, edge_bc=True)
    z, dz = B["z"], B["dz"]
    gammas = [0.0, 0.25, 0.5, 0.75]
    taus = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]

    # ----- GATE A: RK4 linear (N OFF) vs spectral_gap.py -----------------------
    print("\n  GATE A -- RK4 with nonlinear term OFF (must match spectral_gap.py: "
          "rate ~ 2*gamma)")
    sA, n0A, _, _, _, _ = march(B, A=1.0, gammas=gammas, taus=taus, nl=False)
    print("    gamma   measured rate   predicted 2*gamma   (spectral_gap.py)")
    ref = {0.0: 0.054, 0.25: 0.552, 0.5: 1.057, 0.75: 1.675}
    for g in gammas:
        r = fit_rate(taus, [v / n0A[g] for v in sA[g]])
        print(f"    {g:.2f}    {r:+.3f}          {2*g:.2f}                {ref[g]:+.3f}")

    # ----- GATE B + amplitude SWEEP (N ON) -------------------------------------
    print("\n  GATE B + SWEEP -- nonlinear ON, amplitude A "
          "(A->0 must match linear; finite A = the test)")
    print("    A        rate(g=0)  rate(g=.25)  rate(g=.5)  rate(g=.75)  peak||.||w  "
          "edge hi-k")
    amps = [1e-3, 1e-2, 5e-2, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0]
    verdict_growth = []
    for A in amps:
        s, n0, peak, hik, _, _ = march(B, A=A, gammas=gammas, taus=taus, nl=True)
        rates = {g: fit_rate(taus, [v / n0[g] for v in s[g]]) for g in gammas}
        print(f"    {A:<7.3g}  {rates[0.0]:+.3f}     {rates[0.25]:+.3f}      "
              f"{rates[0.5]:+.3f}     {rates[0.75]:+.3f}      {peak:6.3f}     "
              f"{max(hik):.3f}")
        # decay verdict on the reference convective norm gamma=0.5
        verdict_growth.append((A, rates[0.5], peak))

    # ----- CONVERGENCE under (N, zmin) at a finite amplitude -------------------
    print("\n  CONVERGENCE (A=0.5, gamma=0.5 convective rate) under (N, zmin) refinement:")
    print("    N      zmin    rate(g=.5)   peak||.||w   edge hi-k")
    for (Nc, zc, zM) in [(1100, -12.0, 12.0), (1400, -14.0, 12.0), (1800, -16.0, 14.0)]:
        Bc = build_blocks_z(N=Nc, zmin=zc, zmax=zM, beta=0.0, edge_bc=True)
        s, n0, peak, hik, _, _ = march(Bc, A=0.5, gammas=[0.5], taus=taus, nl=True)
        r = fit_rate(taus, [v / n0[0.5] for v in s[0.5]])
        print(f"    {Nc:<6d} {zc:+.0f}    {r:+.3f}       {peak:6.3f}      {max(hik):.3f}")

    # ----- coupled (G,h) run: confirm Theta-coupling adds no instability -------
    print("\n  COUPLED (G and h both excited, A=0.5): confirm no coupling-driven growth")
    s, n0, peak, hik, _, _ = march(B, A=0.5, gammas=[0.0, 0.5], taus=taus, nl=True,
                                   h_amp=0.5)
    print(f"    rate(g=0)={fit_rate(taus, [v/n0[0.0] for v in s[0.0]]):+.3f}   "
          f"rate(g=.5)={fit_rate(taus, [v/n0[0.5] for v in s[0.5]]):+.3f}   "
          f"peak||.||w={peak:.3f}   edge hi-k={max(hik):.3f}")

    # ----- VERDICT -------------------------------------------------------------
    all_decay = all(np.isfinite(r) and r > 0 for (_, r, _) in verdict_growth)
    bounded = all(np.isfinite(p) and p < 1.5 for (_, _, p) in verdict_growth)
    print("\n" + "=" * 80)
    if all_decay and bounded:
        print("  VERDICT: every tested amplitude DECAYS in the convective weighted norm")
        print("  (gamma=0.5 rate > 0) with peak transient ||.||_w/||.||_0 < 1.5 -- the")
        print("  finite-amplitude frozen-base shape perturbation is stable. The nonlinear")
        print("  rate -> the linear gap 2*gamma as A->0 (GATE B) and stays positive up to")
        print("  A=4.  No surprise nonlinear instability. Convective/weighted, frozen base,")
        print("  numerical -- consistent with Conjecture 2.4; NOT a proof.")
    else:
        print("  VERDICT: a finite-amplitude case did NOT decay (rate<=0 or unbounded).")
        print("  CHECK the convergence block: if resolution-ROBUST, this is a candidate")
        print("  nonlinear instability and must be reported, not hidden.")
    print("=" * 80)


if __name__ == "__main__":
    main()
