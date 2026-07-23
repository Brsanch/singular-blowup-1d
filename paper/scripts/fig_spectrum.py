#!/usr/bin/env python3
"""
Figure 1: linear spectrum of M_GG (edge-resolved, inflow BC).
One-command:  python3 paper/scripts/fig_spectrum.py
Output:       paper/figures/fig_spectrum.pdf  (+ .png)
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
from stability_edge_z import build_blocks_z, spectrum  # noqa: E402

OUT = os.path.join(ROOT, "paper", "figures")
os.makedirs(OUT, exist_ok=True)


def main() -> None:
    B = build_blocks_z(N=1000, zmin=-14.0, zmax=8.0, edge_bc=True)
    vGG, _ = spectrum(B["M_GG"])
    vhh, _ = spectrum(B["M_hh"])
    max_re = float(vGG.real.max())
    npos = int((vGG.real > 1e-6).sum())

    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    ax.scatter(vGG.real, vGG.imag, s=10, c="C0", alpha=0.75, label=r"$M_{GG}$", zorder=3)
    ax.scatter(vhh.real, vhh.imag, s=8, c="C1", alpha=0.45, marker="x", label=r"$M_{hh}$")
    ax.axvline(0.0, color="k", lw=0.8, ls="--")
    ax.axhline(0.0, color="0.7", lw=0.5)
    ax.set_xlabel(r"$\mathrm{Re}\,\lambda$")
    ax.set_ylabel(r"$\mathrm{Im}\,\lambda$")
    ax.set_title(
        rf"Edge-resolved spectrum (inflow BC): "
        rf"$\max\mathrm{{Re}}={max_re:+.4f}$, $\#(\mathrm{{Re}}>0)={npos}$"
    )
    # focus on right half-plane neighborhood
    ax.set_xlim(-4.5, 0.4)
    ax.set_ylim(-25, 25)
    ax.legend(loc="upper left", fontsize=9)
    ax.text(
        0.98,
        0.02,
        "KEY maxRe_GG={:.6g} npos={}".format(max_re, npos),
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=7,
        family="monospace",
        color="0.3",
    )
    fig.tight_layout()
    pdf = os.path.join(OUT, "fig_spectrum.pdf")
    png = os.path.join(OUT, "fig_spectrum.png")
    fig.savefig(pdf, dpi=200)
    fig.savefig(png, dpi=160)
    print(f"WROTE {pdf}")
    print(f"WROTE {png}")
    print(f"KEY maxRe_GG={max_re:.8g} npos={npos}")
    # fail-able: paper claims no Re>0
    if npos != 0 or max_re > 0.05:
        print("FAIL: unstable mode present")
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    main()
