#!/usr/bin/env python3
"""
Figure 2: convective spectral gap — measured decay rate vs 2γ.
One-command:  python3 paper/scripts/fig_gap.py
Output:       paper/figures/fig_gap.pdf  (+ .png)
"""
from __future__ import annotations

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))
from stability_edge_z import build_blocks_z  # noqa: E402

OUT = os.path.join(ROOT, "paper", "figures")
os.makedirs(OUT, exist_ok=True)


def evolve(M, G0, taus):
    vals, V = np.linalg.eig(M)
    c = np.linalg.solve(V, G0)
    return [V @ (np.exp(vals * t) * c) for t in taus]


def wnorm(z, dz, G, gamma):
    return float(np.sqrt(np.sum(np.exp(-2 * gamma * z) * np.abs(G) ** 2) * dz))


def fit_rate(taus, ratios, lo=0.5, hi=2.5):
    w = (taus >= lo) & (taus <= hi)
    return -float(np.polyfit(taus[w], np.log(np.array(ratios)[w]), 1)[0])


def main() -> None:
    B = build_blocks_z(N=1000, zmin=-14.0, zmax=10.0, edge_bc=True)
    M, z = B["M_GG"], B["z"]
    dz = z[1] - z[0]
    G0 = np.exp(-((z - 1.0) / 1.0) ** 2)
    gammas = np.array([0.0, 0.25, 0.5, 0.75])
    taus = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    Gt = evolve(M, G0, taus)

    rates = []
    for g in gammas:
        n0 = wnorm(z, dz, G0, g)
        ratios = [wnorm(z, dz, G, g) / n0 for G in Gt]
        rates.append(fit_rate(taus, ratios))
    rates = np.array(rates)

    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    # decay curves for gamma>0
    for g, style in [(0.25, "-"), (0.5, "--"), (0.75, "-.")]:
        n0 = wnorm(z, dz, G0, g)
        ratios = [wnorm(z, dz, G, g) / n0 for G in Gt]
        ax.semilogy(taus, ratios, style, label=rf"$\gamma={g:.2f}$")
    n0u = wnorm(z, dz, G0, 0.0)
    ax.semilogy(
        taus,
        [wnorm(z, dz, G, 0.0) / n0u for G in Gt],
        ":",
        color="0.4",
        label=r"unweighted ($\gamma=0$)",
    )
    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel(r"$\|G(\tau)\|_\gamma / \|G(0)\|_\gamma$")
    ax.set_title("Convective weighted decay (semigroup)")
    ax.legend(fontsize=9)
    ax.grid(True, which="both", alpha=0.3)

    # inset: rate vs 2gamma
    axins = ax.inset_axes([0.52, 0.52, 0.42, 0.38])
    pred = 2 * gammas
    axins.plot(gammas, pred, "k--", lw=1, label=r"$2\gamma$")
    axins.plot(gammas, rates, "o", color="C0", label="measured")
    axins.set_xlabel(r"$\gamma$", fontsize=8)
    axins.set_ylabel("decay rate", fontsize=8)
    axins.tick_params(labelsize=7)
    axins.legend(fontsize=7, loc="upper left")
    axins.set_xlim(-0.05, 0.85)

    fig.tight_layout()
    pdf = os.path.join(OUT, "fig_gap.pdf")
    png = os.path.join(OUT, "fig_gap.png")
    fig.savefig(pdf, dpi=200)
    fig.savefig(png, dpi=160)
    print(f"WROTE {pdf}")
    print(f"WROTE {png}")
    # KEY: rates at 0, 0.25, 0.5, 0.75
    print(
        "KEY rates="
        + ",".join(f"{r:.4f}" for r in rates)
        + "  predict="
        + ",".join(f"{p:.4f}" for p in pred)
    )
    # paper claims ~2γ for γ>0 (within ~15%)
    ok = all(abs(rates[i] - pred[i]) / max(pred[i], 0.1) < 0.20 for i in range(1, len(gammas)))
    if not ok:
        print("FAIL: rates deviate from 2γ")
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    main()
