#!/usr/bin/env python3
"""
Generate publication-quality vector PDF figures for Target B:
  paper/figures_gclm/fig1_gclm_branches_and_spectra.pdf
  paper/figures_gclm/fig2_gclm_coercivity_and_closure.pdf
  paper/figures_gclm/fig3_gclm_dns_selection.pdf

All data are computed directly from the certified src/gclm_singular/ library.
"""

import os
import sys
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))

from gclm_singular.profiles import (  # noqa: E402
    mu_of_a,
    omega_bar,
    advection_bar,
    interior_stagnation_point,
    branch2_omega_bar,
    branch2_mellin_symbol_hl,
    mellin_symbol_c,
)
from gclm_singular.operator import (  # noqa: E402
    build_operator_z,
    build_operator_zeta,
    spectrum,
    _d1_upwind2,
    hilbert_support_nonuniform,
)
from gclm_singular.dns import run_case22_dns  # noqa: E402

OUT_DIR = os.path.join(ROOT, "paper", "figures_gclm")
os.makedirs(OUT_DIR, exist_ok=True)

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.size": 9.5,
        "axes.labelsize": 10,
        "axes.titlesize": 10.5,
        "legend.fontsize": 8.2,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linestyle": "--",
        "lines.linewidth": 1.6,
        "figure.dpi": 150,
    }
)


def make_fig1():
    print("Generating fig1_gclm_branches_and_spectra.pdf ...")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.0))

    # Left panel: Two-branch phase diagram & transcritical exchange at a = -1
    a_line = np.linspace(-2.5, -0.05, 300)
    lam_I = -(1.0 + a_line)
    lam_II = 1.0 + a_line
    ax1.axvspan(-2.5, -1.0, color="#e8f1fa", alpha=0.55, label=r"$a<-1$: Branch I ($L^2_{\mathrm{loc}}$ attractor)")
    ax1.axvspan(-1.0, 0.0, color="#fdf2e9", alpha=0.55, label=r"$-1<a<0$: Branch II ($\beta=\frac{1}{2}$ fixed-cutoff)")
    ax1.axhline(0.0, color="k", lw=0.8, ls="-")
    ax1.axvline(-1.0, color="k", lw=0.9, ls=":")
    ax1.plot(a_line, lam_I, color="#1f5aa6", lw=2.0, label=r"Branch I mode $\lambda_{\mathrm{I},*}(a)=-(1+a)$")
    ax1.plot(a_line, lam_II, color="#c0392b", lw=2.0, ls="--", label=r"Branch II mode $\lambda_{\mathrm{II},*}(a)=1+a$")
    ax1.scatter([-1.0], [0.0], color="k", s=45, zorder=5)
    ax1.annotate(
        r"Transcritical crossing ($a=-1$)" + "\n" + r"$\bar\Omega_{\mathrm{I},-1}=\bar\Omega_{\mathrm{II},-1}$",
        xy=(-1.0, 0.0),
        xytext=(-1.95, 0.72),
        arrowprops=dict(arrowstyle="->", lw=0.9, color="k"),
        fontsize=8.2,
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.7", alpha=0.9),
    )
    ax1.set_xlim(-2.5, 0.0)
    ax1.set_ylim(-1.6, 1.6)
    ax1.set_xlabel(r"Advection parameter $a$")
    ax1.set_ylabel(r"Exact partner eigenvalue $\lambda_*(a)$")
    ax1.set_title(r"(a) Transcritical Exchange of Branches I & II")
    ax1.legend(loc="lower left", framealpha=0.92)

    # Right panel: Computed linear stability spectrum lambda(a) in Gauge (a) and Gauge (b)
    a_list = [-2.5, -2.0, -1.5, -1.25, -1.1, -1.0, -0.9, -0.75, -0.5, -0.25, -0.1]
    lam0_a = []
    lam1_b = []
    for a in a_list:
        if a <= -1.0:
            ops_a = build_operator_z(a, N=350, zmin=-12.0, zmax=10.0, gauge="a", advection_scheme="upwind2")
            ops_b = build_operator_z(a, N=350, zmin=-12.0, zmax=10.0, gauge="b", advection_scheme="upwind2")
        else:
            ops_a = build_operator_zeta(a, N=350, zetamax=12.0, gauge="a", advection_scheme="upwind2")
            ops_b = build_operator_zeta(a, N=350, zetamax=12.0, gauge="b", advection_scheme="upwind2")
        ev_a, _ = spectrum(ops_a["M"])
        ev_b, _ = spectrum(ops_b["M"])
        lam0_a.append(float(np.real(ev_a[0])))
        # First non-zero mode in Gauge (b) (skip near-zero scaling mode)
        nz = [float(np.real(e)) for e in ev_b if abs(e) > 0.05]
        lam1_b.append(nz[0] if nz else float(np.real(ev_b[0])))

    a_arr = np.array(a_list)
    ax2.axhline(0.0, color="k", lw=0.8)
    ax2.axvline(-1.0, color="k", lw=0.9, ls=":")
    a_right = np.linspace(-1.0, -0.05, 100)
    ax2.plot(a_right, -(1.0 + a_right), color="#1f5aa6", lw=1.4, ls="--", label=r"Exact $\lambda_{\mathrm{I},*}(a)=-(1+a)$")
    ax2.plot(a_arr, lam0_a, "o-", color="#1f5aa6", ms=5.0, label=r"Gauge (a) leading $\operatorname{Re}\lambda_0(a)$")
    ax2.plot(a_arr, lam1_b, "s-", color="#d35400", ms=4.8, label=r"Gauge (b) gap $\operatorname{Re}\lambda_1(a)$ (non-zero)")
    ax2.set_xlim(-2.6, 0.0)
    ax2.set_ylim(-1.85, 0.25)
    ax2.set_xlabel(r"Advection parameter $a$")
    ax2.set_ylabel(r"Leading eigenvalue real part $\operatorname{Re}\lambda$")
    ax2.set_title(r"(b) Branch I Linear Spectrum ($\#\operatorname{Re}\lambda>0=0$)")
    ax2.legend(loc="lower right", framealpha=0.92)

    fig.tight_layout()
    out = os.path.join(OUT_DIR, "fig1_gclm_branches_and_spectra.pdf")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out}")


def make_fig2():
    print("Generating fig2_gclm_coercivity_and_closure.pdf ...")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.0))

    # Left panel: Plancherel coercivity symbols (Lemma 2.2 & Theorem 4)
    xi = np.linspace(0.0, 4.0, 600)
    # Branch II at sigma_1 = 5/8 (sigma_- = 3/8), k = 0, 1, 2 for a = -1 and a = -2 (normalized by (1+xi^2)^k)
    r58 = 0.625 + 1j * xi
    for a, col, ls in [(-1.0, "#1f5aa6", "-"), (-2.0, "#c0392b", "--")]:
        tc = branch2_mellin_symbol_hl(r58, a)
        for k, alpha in [(0, 1.0), (2, 0.65)]:
            sym = np.real(((r58 - 0.5) ** k) * (np.conj(r58) ** k) * tc) / ((1.0 + xi ** 2) ** k)
            lbl = rf"Branch II $a={a:g}$, $\sigma_1=\frac{{5}}{{8}}$, $k={k}$"
            ax1.plot(xi, sym, color=col, ls=ls, alpha=alpha, lw=1.6 if k == 0 else 1.2, label=lbl)

    # Branch I Theorem 4 at a = -0.5, mu = 2/3, sigma_1 = mu + 0.25, k = 0, 1
    a_b1 = -0.5
    mu_b1 = mu_of_a(a_b1)
    r_b1 = (mu_b1 + 0.25) + 1j * xi
    ca = mellin_symbol_c(mu_b1 - r_b1, a_b1)
    ax1.plot(xi, np.real(ca), color="#1e8449", lw=1.8, ls="-.", label=r"Branch I Thm 4 $a=-\frac{1}{2}$, $k=0$")
    ax1.plot(
        xi,
        np.real((r_b1 - mu_b1) * np.conj(r_b1) * ca) / (1.0 + xi ** 2),
        color="#1e8449",
        lw=1.3,
        ls=":",
        label=r"Branch I Thm 4 $a=-\frac{1}{2}$, $k=1$",
    )
    ax1.axhline(0.0, color="k", lw=0.8)
    ax1.set_xlim(0.0, 4.0)
    ax1.set_ylim(-2.2, 0.15)
    ax1.set_xlabel(r"Mellin frequency $\xi$")
    ax1.set_ylabel(r"Normalized $\operatorname{Re}[\mathrm{symbol}]/(1+\xi^2)^k$")
    ax1.set_title(r"(a) Unconditional Gauge-(c) Plancherel Coercivity ($<0$)")
    ax1.legend(loc="lower right", fontsize=7.6, framealpha=0.92)

    # Right panel: Nonlinear Gauge-(c) closure on non-analytic C^{3,1} compact bump
    N = 500
    z = np.linspace(-12.0, 14.0, N)
    dz = z[1] - z[0]
    v = np.exp(z)
    Hv = hilbert_support_nonuniform(v, v, dz, v_left=0.0, v_right=v[-1] + 0.5 * (v[-1] - v[-2]), grid_support=True)
    Ht0_op = (v[:, None] * Hv) * (v ** -1.5)[None, :]
    EA = np.zeros((N, N))
    for i in range(1, N):
        tw = np.full(i + 1, dz)
        tw[0] *= 0.5
        tw[i] *= 0.5
        EA[i, : i + 1] = np.exp(z[: i + 1] - z[i]) * tw
    EA[:, 0] += 0.5 * np.exp(z[0] - z)
    Qt0_op = EA @ Ht0_op
    w_z = np.full(N, dz)
    w_z[0] *= 0.5
    w_z[-1] *= 0.5
    ht1_row = -(1.0 / np.pi) * np.exp(-0.5 * z) * w_z
    vg_row = -(1.0 / np.pi) * np.exp(0.5 * z) * np.log1p(np.exp(-z)) * w_z
    Dz = _d1_upwind2(N, dz, np.full(N, 2.0))
    w14 = np.exp(-0.5 * z)
    w38 = np.exp(-0.75 * z)
    w58 = np.exp(-1.25 * z)
    m_int = z >= -10.0

    ax_ins = ax2.inset_axes([0.48, 0.50, 0.48, 0.44])
    colors = {-0.5: "#8e44ad", -1.0: "#1f5aa6", -1.5: "#d35400", -2.0: "#c0392b"}
    for a, dt in [(-0.5, 0.001), (-1.0, 0.002), (-1.5, 0.004), (-2.0, 0.004)]:
        G = 0.2 * np.maximum(0.0, 1.0 - (0.5 * z) ** 2) ** 4
        G[0] = 0.0
        wm = w38 if a < -1.6 else w14

        def rhs(g):
            Ht0 = Ht0_op @ g
            Qt0 = Qt0_op @ g
            ht1 = float(ht1_row @ g)
            Vg = float(vg_row @ g)
            r = -(2.0 + Qt0 + ht1 - Vg) * (Dz @ g) + ((1.0 / a) * Ht0 + 0.5 * Qt0) * (1.0 + g)
            r[0] = 0.0
            cl = 2.0 - Vg
            cw = -0.5 * cl - (1.0 / a + 0.5) * ht1
            return r, cl, cw

        T_end = 4.0
        nsteps = int(round(T_end / dt))
        samp = max(1, int(round(0.1 / dt)))
        ts, Es, perrs = [], [], []
        for st in range(nsteps + 1):
            if st % samp == 0:
                f1 = np.gradient(G, dz)
                f2 = np.gradient(f1, dz)
                f3 = np.gradient(f2, dz)
                E_val = float(
                    np.sum(wm[m_int] * (G[m_int] ** 2 + f1[m_int] ** 2 + f2[m_int] ** 2)) * dz
                    + np.sum(w58[m_int] * (G[m_int] ** 2 + f1[m_int] ** 2 + f2[m_int] ** 2 + f3[m_int] ** 2)) * dz
                )
                _, cl, cw = rhs(G)
                ts.append(st * dt)
                Es.append(E_val)
                perrs.append(abs(cl - 2.0) + abs(cw + 1.0))
            if st < nsteps:
                k1, *_ = rhs(G)
                k2, *_ = rhs(G + 0.5 * dt * k1)
                k3, *_ = rhs(G + 0.5 * dt * k2)
                k4, *_ = rhs(G + dt * k3)
                G = G + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
                G[0] = 0.0
        ax2.semilogy(ts, Es, color=colors[a], lw=1.7, label=rf"$a={a:g}$ (Branch I at $a=-1$)" if a == -1.0 else rf"$a={a:g}$")
        ax_ins.semilogy(ts, perrs, color=colors[a], lw=1.2)

    t_ref = np.linspace(0.0, 4.0, 50)
    ax2.semilogy(t_ref, 3.0 * np.exp(-1.5 * t_ref), "k--", lw=1.0, alpha=0.6, label=r"ref. $e^{-1.5\tau}$")
    ax2.set_xlabel(r"Rescaled time $\tau$")
    ax2.set_ylabel(r"Sobolev energy $\mathcal{E}(\tau)$")
    ax2.set_title(r"(b) Nonlinear Gauge-(c) Closure (Non-Analytic $C^{3,1}$ Data)")
    ax2.legend(loc="lower left", fontsize=7.8, framealpha=0.92)
    ax_ins.set_title(r"$|c_l-2|+|c_\omega+1|$", fontsize=7.5, pad=2)
    ax_ins.tick_params(labelsize=6.8)
    ax_ins.set_xlim(0.0, 4.0)

    fig.tight_layout()
    out = os.path.join(OUT_DIR, "fig2_gclm_coercivity_and_closure.pdf")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out}")


def make_fig3():
    print("Generating fig3_gclm_dns_selection.pdf ...")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.0))

    # Left panel: Universal convergence at a = -2.0 to Branch I vs Branch II
    res_exp = run_case22_dns(-2.0, init_type="exp_inv", N=450, c_stretch=1e-4, rhomax=29.0, T_end=25.0)
    res_b2 = run_case22_dns(-2.0, init_type="branch2", N=450, c_stretch=1e-4, rhomax=29.0, T_end=25.0)
    v = res_exp["v"]
    X = res_exp["X"]
    m = (v >= 5e-3) & (v <= 5e2)
    Om_I = np.abs(omega_bar(X[m], -2.0))
    Om_II = np.abs(branch2_omega_bar(X[m], -2.0, unit_H0=True))

    ax1.loglog(v[m], Om_I, "k-", lw=2.2, label=r"Exact Branch I $|\bar\Omega_{\mathrm{I},-2}|$ ($\beta=\frac{1}{3}$)")
    ax1.loglog(v[m], Om_II, color="#c0392b", lw=1.6, ls="--", label=r"Exact Branch II $|\bar\Omega_{\mathrm{II},-2}|$ ($\beta=\frac{1}{2}$)")
    ax1.loglog(v[m], np.abs(res_exp["Om"][m]), color="#1f5aa6", lw=1.8, ls=":", label=r"DNS from smooth $e^{-1/v}$")
    ax1.loglog(v[m], np.abs(res_b2["Om"][m]), color="#1e8449", lw=1.4, ls="-.", label=r"DNS from reg. Branch II IC")
    ax1.axvline(0.125, color="0.4", lw=0.9, ls=":")
    ax1.annotate(
        r"$v_*=X_*-1=\frac{1}{8}$",
        xy=(0.125, 1.5),
        xytext=(0.35, 3.0),
        arrowprops=dict(arrowstyle="->", lw=0.8, color="0.3"),
        fontsize=8.0,
    )
    ax1.set_xlabel(r"Distance from edge $v = X - 1$")
    ax1.set_ylabel(r"Profile magnitude $|\Omega(1+v)|$")
    ax1.set_title(r"(a) Universal Convergence to Branch I ($a=-2.0$)")
    ax1.legend(loc="lower left", fontsize=7.6, framealpha=0.92)

    # Right panel: Effective scaling exponent gamma(a) & singularity slope beta(a) selected by Case 2.2 DNS
    res_15 = run_case22_dns(-1.5, init_type="exp_inv", N=450, c_stretch=1e-4, rhomax=29.0, T_end=25.0)
    res_25 = run_case22_dns(-2.5, init_type="exp_inv", N=450, c_stretch=1e-4, rhomax=29.0, T_end=25.0)
    res_05 = run_case22_dns(-0.5, init_type="exp_inv", N=450, c_stretch=1e-4, rhomax=29.0, T_end=15.0)

    a_dense = np.linspace(-2.6, -0.05, 200)
    ax2.plot(a_dense, 1.0 - a_dense, color="#1f5aa6", lw=2.0, label=r"Branch I exact $\gamma_{\mathrm{I}}(a)=1-a$")
    ax2.plot(a_dense, np.full_like(a_dense, 2.0), color="#c0392b", lw=1.8, ls="--", label=r"Branch II exact $\gamma_{\mathrm{II}}(a)=2$")
    ax2.axvline(-1.0, color="k", lw=0.9, ls=":")
    ax2.axvspan(-2.6, -1.0, color="#e8f1fa", alpha=0.45)
    ax2.axvspan(-1.0, 0.0, color="#fdf2e9", alpha=0.45)

    dns_pts = [(-2.5, res_25["gamma"]), (-2.0, res_exp["gamma"]), (-1.5, res_15["gamma"]), (-0.5, res_05["gamma"])]
    ax2.scatter(
        [p[0] for p in dns_pts[:3]],
        [p[1] for p in dns_pts[:3]],
        color="#1f5aa6",
        edgecolor="white",
        s=65,
        zorder=6,
        label=r"Case 2.2 DNS ($a<-1 \to$ Branch I)",
    )
    ax2.scatter(
        [dns_pts[3][0]],
        [dns_pts[3][1]],
        color="#c0392b",
        marker="D",
        edgecolor="white",
        s=65,
        zorder=6,
        label=r"Fixed-cutoff DNS ($a=-0.5 \to$ Branch II)",
    )
    ax2.annotate(
        rf"$\gamma={res_05['gamma']:.2f},\ \beta={res_05['beta_fit']:.2f}$",
        xy=(-0.5, res_05["gamma"]),
        xytext=(-0.92, 2.55),
        arrowprops=dict(arrowstyle="->", lw=0.8, color="#c0392b"),
        fontsize=8.0,
        color="#c0392b",
    )
    ax2.set_xlim(-2.6, -0.05)
    ax2.set_ylim(1.0, 3.85)
    ax2.set_xlabel(r"Advection parameter $a$")
    ax2.set_ylabel(r"Selected shrinking exponent $\gamma = -c_l/c_\omega$")
    ax2.set_title(r"(b) $L^2$-Criticality Selection Across $a=-1$")
    ax2.legend(loc="upper right", fontsize=7.5, framealpha=0.92)

    fig.tight_layout()
    out = os.path.join(OUT_DIR, "fig3_gclm_dns_selection.pdf")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out}")


if __name__ == "__main__":
    make_fig1()
    make_fig2()
    make_fig3()
    print("All 3 Target B figures generated successfully.")
