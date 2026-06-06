#!/usr/bin/env python3
"""
STRESS TEST of the self-consistent co-evolving attractor (falsification-first).

coevolve_projected.py showed ONE perturbation (A=0.15 bulk bump) relaxes to
(Omega_bar, c_l=2, c_omega=-1). Claim under test: the singular profile is a nonlinear
attractor. Try to BREAK it: sweep amplitude and perturbation SHAPE (bulk / near-edge /
far-field / high-k oscillatory / Theta-dominant) in the fully self-consistent Theta-
balance march, and check each relaxes (||g||_w decays, c_l->2, c_omega->-1) rather than
blowing up or drifting. A resolution-robust failure would be a real counterexample.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from coevolve_projected import rk4_step, gauge_of  # noqa: E402
from coevolve_moving_frame import build, wnorm, fit_rate  # noqa: E402


def march_ic(B, G0, h0, T=4.0, taus=(0.0, 0.5, 1.0, 2.0, 3.0, 4.0)):
    G = G0.copy(); h = h0.copy(); G[0] = 0.0; h[0] = 0.0
    xi = 0.0
    dt = 0.25 * B["dz"]
    _, c_l, c_om = gauge_of(B, G, h, "theta_balance")
    series = [(0.0, c_l, c_om, wnorm(B, G), float(np.max(np.abs(G))))]
    tau, ti = 0.0, 1
    targets = list(taus)
    while ti < len(targets):
        while tau < targets[ti] - 1e-12:
            G, h, xi, c_l, c_om = rk4_step(B, G, h, xi, min(dt, targets[ti] - tau),
                                           mode="theta_balance")
            tau += min(dt, targets[ti] - tau)
            if not np.all(np.isfinite(G)):
                series.append((tau, np.nan, np.nan, np.inf, np.inf)); return series
        series.append((targets[ti], c_l, c_om, wnorm(B, G), float(np.max(np.abs(G)))))
        ti += 1
    return series


def shapes(z):
    return {
        "bulk bump (z=1)":     np.exp(-((z - 1.0) / 1.0) ** 2),
        "near-edge (z=-4)":    np.exp(-((z + 4.0) / 1.0) ** 2),
        "far-field (z=6)":     np.exp(-((z - 6.0) / 1.0) ** 2),
        "high-k osc (z=1)":    np.exp(-((z - 1.0) / 1.5) ** 2) * np.cos(6.0 * z),
        "two-bump":            np.exp(-((z - 0.0) / 0.8) ** 2) - np.exp(-((z - 3.0) / 0.8) ** 2),
    }


def verdict(series):
    last = series[-1]
    if not np.isfinite(last[3]):
        return "BLOWS UP"
    rate = fit_rate([s[0] for s in series], [s[3] for s in series])
    return (f"relax: ||g||_w {series[0][3]:.2e}->{last[3]:.2e}  rate {rate:+.2f}  "
            f"c_l={last[1]:+.5f} c_om={last[2]:+.5f}")


def main():
    print("=" * 84)
    print("  STRESS TEST -- self-consistent co-evolving attractor (try to break relaxation)")
    print("=" * 84)
    B = build(N=1400, zmin=-12.0, zmax=12.0)
    z = B["z"]

    print("\n  AMPLITUDE sweep (bulk bump on Omega-shape):")
    print("    A        verdict")
    for A in [0.05, 0.15, 0.3, 0.6, 1.0, 2.0]:
        s = march_ic(B, A * np.exp(-((z - 1.0) / 1.0) ** 2), np.zeros(B["n"]))
        print(f"    {A:<7.3g}  {verdict(s)}")

    print("\n  SHAPE sweep (A=0.2 on Omega-shape):")
    print("    shape                  verdict")
    for name, sh in shapes(z).items():
        s = march_ic(B, 0.2 * sh, np.zeros(B["n"]))
        print(f"    {name:22s} {verdict(s)}")

    print("\n  THETA-dominant perturbation (Omega bump 0.05 + big Theta bump 0.5):")
    s = march_ic(B, 0.05 * np.exp(-((z - 1.0) / 1.0) ** 2),
                 0.5 * np.exp(-((z - 1.5) / 1.0) ** 2))
    print(f"    {verdict(s)}")

    print("\n  CONVERGENCE of the most aggressive case (A=1.0 bulk) under (N,zmin):")
    print("    N      zmin   verdict")
    for (Nc, zc_) in [(1000, -10.0), (1400, -12.0), (1800, -14.0)]:
        Bc = build(N=Nc, zmin=zc_, zmax=12.0)
        s = march_ic(Bc, 1.0 * np.exp(-((Bc["z"] - 1.0) / 1.0) ** 2), np.zeros(Bc["n"]))
        print(f"    {Nc:<6d} {zc_:+.0f}   {verdict(s)}")

    print("\n" + "=" * 84)
    print("  READ: every amplitude and shape that decays (||g||_w down, c_l->2, c_omega->-1,")
    print("  no blowup) strengthens 'the singular profile is a nonlinear attractor of the")
    print("  self-consistent co-evolving flow'. A RESOLUTION-ROBUST blowup/drift would be a")
    print("  genuine counterexample (report it, do not hide). Convective/weighted, numerical.")
    print("=" * 84)


if __name__ == "__main__":
    main()
