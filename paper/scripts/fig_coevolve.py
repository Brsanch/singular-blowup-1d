#!/usr/bin/env python3
"""
Figure 3: co-evolving march — (c_l, c_ω) → (2, −1).
One-command:  python3 paper/scripts/fig_coevolve.py
Output:       paper/figures/fig_coevolve.pdf  (+ .png)
"""
from __future__ import annotations

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "validate"))
sys.path.insert(0, os.path.join(ROOT, "src"))
from coevolve_moving_frame import build  # noqa: E402
from coevolve_projected import gauge_of, rk4_step, wnorm  # noqa: E402

OUT = os.path.join(ROOT, "paper", "figures")
os.makedirs(OUT, exist_ok=True)


def march_traj(B, A=0.15, T=4.0, n_out=40):
    z = B["z"]
    G = A * np.exp(-((z - 1.0) / 1.0) ** 2)
    G[0] = 0.0
    h = np.zeros_like(G)
    xi = 0.0
    dt = 0.25 * B["dz"]
    taus = np.linspace(0.0, T, n_out)
    rows = []
    _, c_l, c_om = gauge_of(B, G, h, "theta_balance")
    rows.append((0.0, c_l, c_om, wnorm(B, G)))
    tau, ti = 0.0, 1
    while ti < len(taus):
        target = taus[ti]
        while tau < target - 1e-12:
            step = min(dt, target - tau)
            G, h, xi, c_l, c_om = rk4_step(B, G, h, xi, step, mode="theta_balance")
            tau += step
        rows.append((tau, c_l, c_om, wnorm(B, G)))
        ti += 1
    return np.array(rows)


def main() -> None:
    # Moderate resolution for figure (gate uses multiple grids; paper cites full gate)
    B = build(N=500, zmin=-12.0, zmax=8.0)
    traj = march_traj(B, A=0.15, T=4.0, n_out=33)
    tau, cl, com, wg = traj.T
    cl_final, com_final = float(cl[-1]), float(com[-1])

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4))

    ax = axes[0]
    ax.plot(tau, cl, "C0-", lw=1.8, label=r"$c_l(\tau)$")
    ax.axhline(2.0, color="C0", ls="--", lw=0.9, alpha=0.7)
    ax.plot(tau, com, "C1-", lw=1.8, label=r"$c_\omega(\tau)$")
    ax.axhline(-1.0, color="C1", ls="--", lw=0.9, alpha=0.7)
    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel("exponents")
    ax.set_title(r"Co-evolving exponents $\to(2,-1)$")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 4)

    ax = axes[1]
    ax.semilogy(tau, wg / wg[0], "C2-", lw=1.8)
    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel(r"$\|G\|_w / \|G\|_w(0)$")
    ax.set_title("Weighted shape residual")
    ax.grid(True, which="both", alpha=0.3)
    ax.set_xlim(0, 4)

    fig.tight_layout()
    pdf = os.path.join(OUT, "fig_coevolve.pdf")
    png = os.path.join(OUT, "fig_coevolve.png")
    fig.savefig(pdf, dpi=200)
    fig.savefig(png, dpi=160)
    print(f"WROTE {pdf}")
    print(f"WROTE {png}")
    print(f"KEY c_l_final={cl_final:.8g} c_omega_final={com_final:.8g}")
    # paper: converge to (2,-1) — figure-resolution: within 1e-2
    if abs(cl_final - 2.0) > 1e-2 or abs(com_final + 1.0) > 1e-2:
        print("FAIL: exponents did not relax to (2,-1)")
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    main()
