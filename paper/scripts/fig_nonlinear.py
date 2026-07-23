#!/usr/bin/env python3
"""
Figure 4: nonlinear finite-amplitude decay (real march from nonlinear_shape_march)
           + coupled resolvent channels.
One-command:  python3 paper/scripts/fig_nonlinear.py
Output:       paper/figures/fig_nonlinear.pdf  (+ .png)
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
sys.path.insert(0, os.path.join(ROOT, "validate"))
from stability_edge_z import build_blocks_z  # noqa: E402
from nonlinear_shape_march import march, fit_rate  # noqa: E402

OUT = os.path.join(ROOT, "paper", "figures")
os.makedirs(OUT, exist_ok=True)


def coupled_gains(B, sigma=0.4, lam=0.3):
    z, dz = B["z"], B["dz"]
    Mgg, Mhh, Mgh = B["M_GG"], B["M_hh"], B["M_Gh"]
    n = len(z)
    I = np.eye(n)

    def wn(F):
        return float(np.sqrt(np.sum(np.exp(-2 * sigma * z) * np.abs(F) ** 2) * dz))

    z0s = np.array([-2.0, -4.0, -6.0, -8.0, -10.0])
    gO, gTh = [], []
    for z0 in z0s:
        Fb = np.exp(-np.cosh((z - z0)) + 1.0)
        nF = wn(Fb)
        G = np.linalg.solve(lam * I - Mgg, Fb)
        gO.append(wn(G) / nF)
        h = np.linalg.solve(lam * I - Mhh, Fb)
        G2 = np.linalg.solve(lam * I - Mgg, Mgh @ h)
        gTh.append(wn(G2) / nF)
    return z0s, np.array(gO), np.array(gTh)


def main() -> None:
    # Match gate resolution for KEY numbers
    B = build_blocks_z(N=1100, zmin=-12.0, zmax=10.0, edge_bc=True)
    gammas = [0.5]  # convective rate reported in paper
    taus = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
    As = np.array([0.1, 0.5, 1.0, 2.0, 4.0])
    rates, peaks = [], []
    for A in As:
        s, n0, peak, hik, _, _ = march(B, A=float(A), gammas=gammas, taus=taus, nl=True)
        r = fit_rate(taus, [v / n0[0.5] for v in s[0.5]])
        # peak ratio at sample times for gamma=0.5
        ratios = [v / n0[0.5] for v in s[0.5]]
        p = float(max(ratios))
        rates.append(r)
        peaks.append(p)
        print(f"  A={A:.1f}: rate={r:.4f}  peak={p:.4f}  hik={hik[-1]:.3f}")
    rates = np.array(rates)
    peaks = np.array(peaks)
    z0s, gO, gTh = coupled_gains(B)

    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.5))
    ax = axes[0]
    ax.plot(As, rates, "o-", color="C0", lw=1.6, label="measured rate")
    ax.axhline(1.0, color="k", ls="--", lw=0.9, label=r"linear gap $2\gamma=1$")
    ax.set_xlabel(r"amplitude $A$")
    ax.set_ylabel(r"weighted decay rate ($\gamma=0.5$)")
    ax.set_title("Nonlinear frozen-base decay")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0.0, max(1.4, float(rates.max()) * 1.15))
    ax.text(
        0.98,
        0.05,
        f"max peak ratio = {peaks.max():.3f}",
        transform=ax.transAxes,
        ha="right",
        fontsize=8,
        color="0.3",
    )

    ax = axes[1]
    ax.semilogy(z0s, gO, "o-", label=r"pure $\Omega$ gain")
    ax.semilogy(z0s, gTh, "s-", label=r"$\Theta\to G$ gain")
    ax.set_xlabel(r"force depth $z_0$")
    ax.set_ylabel(r"weighted resolvent gain")
    ax.set_title(r"Coupled channels ($\lambda=0.3$, $\sigma=0.4$)")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    pdf = os.path.join(OUT, "fig_nonlinear.pdf")
    png = os.path.join(OUT, "fig_nonlinear.png")
    fig.savefig(pdf, dpi=200)
    fig.savefig(png, dpi=160)
    print(f"WROTE {pdf}")
    print(f"WROTE {png}")
    print(
        "KEY rates="
        + ",".join(f"{r:.4f}" for r in rates)
        + " peaks="
        + ",".join(f"{p:.4f}" for p in peaks)
        + f" gTh_deep={gTh[-1]:.4f} gO_deep={gO[-1]:.4f}"
    )
    # Full gate: peak=1.000 at sample times; allow 1.02 for moderate N
    if peaks.max() > 1.02:
        print("FAIL: peak transient growth > 1.02")
        sys.exit(1)
    if np.any(rates < 0.5):
        print("FAIL: rate far below linear gap")
        sys.exit(1)
    if gTh[-1] > 3.0:
        print("FAIL: Theta->G blowup")
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    main()
