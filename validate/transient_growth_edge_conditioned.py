#!/usr/bin/env python3
"""
SETTLE transient growth with the EDGE-CONDITIONED operator.

findings_transient_growth.md left one caveat open: PHYSICAL perturbations have G=1 (no
transient growth), but the WORST-CASE operator norm sup_tau ||e^{M tau}|| crept up with edge
DEPTH (3.08 -> 3.91) -- a grid-scale (Y-1)^{-1/2} edge spike -- and "a function-analytic
statement would need the precise weighted space and an edge cutoff." We now have both
(energy_estimate_search.py): the bounded-at-edge anisotropic weight + the edge cutoff where
the bulk energy abscissa is negative.

KEY FACT: in a norm with weight W, the numerical abscissa (log-norm)
    mu_W(M) = lambda_max( sym_W(M) ) = lambda_max( (M~ + M~^T)/2 ),  M~ = sqrt(W) M sqrt(W)^{-1}
bounds transient growth EXACTLY:  ||e^{M tau}||_W <= e^{mu_W tau}.  So mu_W <= 0  =>
||e^{M tau}||_W <= 1 for all tau  =>  NO transient growth (monotone contraction) in <.>_W.

On z > z_cut with the anisotropic weight W, we compare:
  (1) the continuum-IBP cutoff-model surrogate mu_an = lambda_max(-diag(phi') + sym_W(1/2 Q_bulk)),
      which replaces the discrete transport symmetric part by the analytic IBP drift -phi'(z)
      for functions vanishing at z_cut and is negative and resolution-independent (-0.101 at
      z_cut=-8, -0.407 at z_cut=-6);
  (2) the actual discrete principal submatrix M_b = -2 D_{z,bulk} + 1/2 Q_bulk, whose unclosed
      centered-difference stencil at z_cut gives mu_disc = lambda_max(sym_W(M_b)) > 0 (+33.4 at
      N=1000) and G_mat = sup_tau ||exp(M_b tau)||_W ~ 112 (so M_b is NOT a discrete matrix
      contraction), while smooth interior Gaussian bumps have G(bulk bump) = 1.0000.
"""
import os
import sys

import numpy as np
from scipy.linalg import expm

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from stability_edge_z import build_blocks_z  # noqa: E402

GAMMA, Z0, S = 0.4, -6.0, 2.0          # the bounded-at-edge anisotropic weight (energy_estimate_search)
TAUS = np.linspace(0.0, 6.0, 31)


def build_bounded(N, zmin, zmax):
    B = build_blocks_z(N=N, zmin=zmin, zmax=zmax, beta=0.0, edge_bc=True)
    M = -2.0 * B["Dz"] + 0.5 * B["Q_op"]   # bounded-forcing operator (= transient_growth.py form)
    return M, B["Dz"], B["Q_op"], B["z"], B["dz"]


def worstcase_standard(M):
    """sup_tau ||e^{M tau}||_2 in the plain norm (the OLD measure: edge-artifact creep)."""
    return max(float(np.linalg.norm(expm(M * t), 2)) for t in TAUS)


def edge_conditioned(Dz, Qop, z, dz, z_cut):
    """Bulk (z>z_cut) operator in the bounded-at-edge anisotropic W-norm.
    The numerical abscissa (log-norm) uses the ANALYTIC transport symmetric part -phi'(z)
    (exact IBP), NOT sym(discrete -2 d_z) which carries a boundary-discretization artifact
    (cap_scaffold). So  mu_an = lambda_max( -diag(phi') + sym_W( (1/2) Q_op ) ).
    Also returns the physical bulk-bump growth via the discrete semigroup (smooth perts do
    not excite the grid artifact -> direct confirmation)."""
    keep = z > z_cut
    zk = z[keep]
    pp = 2.0 * GAMMA * 0.5 * (1.0 + np.tanh((zk - Z0) / S))   # phi'(z): 0 at edge -> 2 gamma bulk
    sw = np.exp(-0.5 * (np.cumsum(pp) * dz))                  # sqrt(W)
    Fb = 0.5 * Qop[np.ix_(keep, keep)]                       # bulk nonlocal forcing
    Ft = (sw[:, None] * Fb) / sw[None, :]                     # W-similarity of the forcing
    mu_an = float(np.linalg.eigvalsh(-np.diag(pp) + 0.5 * (Ft + Ft.T))[-1])
    # physical bulk-bump transient growth and actual discrete submatrix norm in W-norm
    Mb = (-2.0 * Dz[np.ix_(keep, keep)] + Fb * 2.0 * 0.5)     # discrete bulk submatrix (-2 d_z + 1/2 Q)
    Mt = (sw[:, None] * Mb) / sw[None, :]
    mu_disc = float(np.linalg.eigvalsh(0.5 * (Mt + Mt.T))[-1])
    Gmat = max(float(np.linalg.norm(expm(Mt * t), 2)) for t in np.r_[0.001, TAUS])
    gb = np.exp(-((zk - 1.5) / 1.2) ** 2); gb = gb / np.linalg.norm(gb)
    Gbump = max(float(np.linalg.norm(expm(Mt * t) @ gb)) for t in TAUS)
    return mu_an, mu_disc, Gmat, Gbump


def main():
    print("=" * 84)
    print("  SETTLE transient growth -- EDGE-CONDITIONED operator (bulk + anisotropic W-norm)")
    print("=" * 84)

    print("\n  UN-conditioned baseline (full operator, PLAIN norm) -- the old edge creep:")
    print("    N      zmin   worst-case sup_tau ||e^{M t}||_2")
    for (N, zmin, zmax) in ((800, -10, 8), (1000, -12, 8), (1400, -12, 8), (1200, -14, 8)):
        M, Dz, Qop, z, dz = build_bounded(N, zmin, zmax)
        print(f"    {N:<6d} {zmin:+.0f}   {worstcase_standard(M):6.2f}")

    print("\n  EDGE-CONDITIONED (bulk z>z_cut, anisotropic W-norm; continuum-IBP surrogate mu_an vs discrete submatrix):")
    print("  mu_an = lambda_max(-diag(phi') + sym_W(1/2 Q_bulk)); G_mat = sup_tau ||exp(M_b tau)||_W")
    print("    N      zmin   z_cut   mu_an (IBP surr)   mu_disc(sym_W)   G_mat(disc)   G(bulk bump)")
    for (N, zmin, zmax) in ((1000, -12, 8), (1400, -12, 8), (2000, -12, 8), (1600, -14, 8)):
        M, Dz, Qop, z, dz = build_bounded(N, zmin, zmax)
        for z_cut in (-8.0, -6.0):
            mu, mu_d, Gm, Gb = edge_conditioned(Dz, Qop, z, dz, z_cut)
            print(f"    {N:<6d} {zmin:+.0f}   {z_cut:+.0f}    {mu:+.5f}          {mu_d:+8.2f}        {Gm:6.4f}        {Gb:6.4f}")

    print("\n" + "=" * 84)
    print("  READ: UN-conditioned plain-norm worst-case creeps with edge depth (2.33 -> 3.91).")
    print("  On z > z_cut, the continuum-IBP cutoff-model surrogate mu_an is negative and")
    print("  resolution-independent (-0.101 at z_cut=-8, -0.407 at -6), and smooth interior")
    print("  bumps satisfy G(bulk bump) = 1.0000. However, the discrete principal submatrix")
    print("  M_b = -2 D_{z,bulk} + 1/2 Q_bulk has an unclosed centered-difference stencil at")
    print("  z_cut, giving mu_disc(sym_W) > 0 (+33.4 to +67.3) and sampled matrix amplification")
    print("  G_mat(disc) ~ 112-113 (1.033 at tau=1e-3 on the top symmetric eigenvector), so")
    print("  mu_an is a continuum cutoff surrogate and discrete matrix contraction is not claimed.")
    print("=" * 84)


if __name__ == "__main__":
    main()
