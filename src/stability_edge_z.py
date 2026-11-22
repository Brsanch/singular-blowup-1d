#!/usr/bin/env python3
"""
Desingularized stability operator in the EDGE-RESOLVING coordinate z = ln(Y-1).

Same operator as src/stability_edge.py (Conjecture 2.4, frozen mu=nu=0), but on a
uniform z-grid (Y = 1 + e^z).  Two payoffs over the equispaced build:
  * z clusters at the edge (Y->1 is z->-inf) and compresses the 1/sqrt(Y) tail, so
    the singular profile is resolved -> no spurious near-edge eigenvalue cloud.
  * transport -abar d_Y = -2(Y-1) d_Y = -2 d_z is CONSTANT-COEFFICIENT, so the
    essential spectrum is a clean vertical line and POINT modes show as isolated
    eigenvalues off it.  A far-field weight w=e^{beta z} (beta<0) pushes the
    essential spectrum strictly left to reveal a spectral GAP (point spectrum is
    weight-invariant).

Weighted variable G = (Y-1)^{1/2} g (g = (Y-1)^{-1/2} G = e^{-z/2} G), so G is
bounded at the edge.  Hilbert transform of the support-only perturbation is done
with a non-equispaced p.v. quadrature (subtraction form, exact on constants),
validated against the closed form H((Y-1)^{-1/2}) = 0 on Y>1.

Block-triangular as before: full spectrum = eig(M_GG) U eig(M_hh) U {c_l-1}.
"""

import numpy as np

C_L = 2.0
C_OMEGA = -1.0


def _dz_central(N, dz):
    """2nd-order central d_z, non-periodic (one-sided at the two ends)."""
    D = np.zeros((N, N))
    for i in range(1, N - 1):
        D[i, i + 1] = 1.0 / (2 * dz)
        D[i, i - 1] = -1.0 / (2 * dz)
    D[0, 0] = -3 / (2 * dz); D[0, 1] = 4 / (2 * dz); D[0, 2] = -1 / (2 * dz)
    D[-1, -1] = 3 / (2 * dz); D[-1, -2] = -4 / (2 * dz); D[-1, -3] = 1 / (2 * dz)
    return D


def _cumz(N, dz):
    """Cumulative trapezoid in z from z_min: (Cz @ f)_i ~ int_{z_min}^{z_i} f dz'."""
    C = np.zeros((N, N))
    for i in range(1, N):
        C[i] = C[i - 1]
        C[i, i] += 0.5 * dz
        C[i, i - 1] += 0.5 * dz
    return C


def hilbert_support_z(Y, Ymax):
    """Non-equispaced line-Hilbert for g supported on Y>1, sampled on a clustered
    Y-grid.  H(g)(Y_i) = (1/pi) pv int_1^{Ymax} g/(Y_i-Y') dY'.  Subtraction form:
    off-diagonal (1/pi) wY_j/(Y_i-Y_j); diagonal absorbs the exact pv of a constant
    so that H is exact on constants (and accurate on smooth g)."""
    n = len(Y)
    wY = np.empty(n)
    wY[1:-1] = 0.5 * (Y[2:] - Y[:-2])
    wY[0] = 0.5 * (Y[1] - Y[0])
    wY[-1] = 0.5 * (Y[-1] - Y[-2])
    d = Y[:, None] - Y[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        K = wY[None, :] / d
    K[np.diag_indices(n)] = 0.0
    H = K / np.pi
    pv = np.log((Y - 1.0) / (Ymax - Y))
    H[np.diag_indices(n)] = (pv - K.sum(axis=1)) / np.pi
    return H


def edge_tail_row(z):
    """Row t with  t @ G = int_0^{v_min} H[g](1+y) dy  (v_min = e^{z_min}), EXACT for
    g = v^{-1/2} G supported on v >= v_min (which is what the grid represents):

        int_0^{v_min} dy / (v' - y) = -ln(1 - v_min/v')   =>
        t @ G = -(1/pi) int_{z_min}^{z_max} e^{z'/2} (-ln(1 - e^{z_min - z'})) G(z') dz'.

    For z' >> z_min the kernel tends to -(1/pi) e^{z_min} e^{-z'/2}, i.e. e^{z_min} times
    the edge functional a = H[g](1).  The log singularity at z' = z_min is integrable;
    weights are exact integrals of the kernel against the piecewise-linear hat functions
    (Gauss-Legendre per cell, adaptive quadrature on the first, singular, cell)."""
    from scipy.integrate import quad
    z = np.asarray(z, dtype=float)
    n = len(z)
    z0 = z[0]

    def k(zp):
        return np.exp(0.5 * zp) * (-np.log1p(-np.exp(z0 - zp)))

    xg, wg = np.polynomial.legendre.leggauss(12)
    W = np.zeros(n)
    a, b = z[1:-1], z[2:]                       # cells j = 1 .. n-2 (smooth)
    half = 0.5 * (b - a)
    zq = a[:, None] + half[:, None] * (xg[None, :] + 1.0)
    kq = k(zq) * wg[None, :] * half[:, None]
    lam = (zq - a[:, None]) / (b - a)[:, None]  # 0 at left node, 1 at right node
    W[1:-1] += (kq * (1.0 - lam)).sum(axis=1)
    W[2:] += (kq * lam).sum(axis=1)
    h = z[1] - z[0]                             # first cell: log singularity at z0
    W[0] += quad(lambda zp: k(zp) * (z[1] - zp) / h, z0, z[1], limit=200)[0]
    W[1] += quad(lambda zp: k(zp) * (zp - z0) / h, z0, z[1], limit=200)[0]
    return -W / np.pi


def build_blocks_z(N=1400, zmin=-14.0, zmax=7.0, beta=0.0, edge_bc=True, edge_tail=True):
    """Assemble the frozen operator on z in [zmin, zmax].

    beta : far-field weight exponent (w = e^{beta z}); beta<0 pushes the essential
           (transport) spectrum strictly left -> reveals the spectral gap.
    edge_bc : impose the inflow/regularity condition G(z_min)=0 (and h(z_min)=0).
           The edge z->-inf is an INFLOW boundary (characteristics move +z); the
           weighted space excludes the homogeneous branch G~(Y-1)^{-lambda/2} that
           blows up at the edge for Re lambda>0.  Without this BC the discretization
           re-admits that inadmissible branch as a spurious z_min-dependent edge
           mode (e.g. Re~+1.7).  Implemented by dropping the z_min DOF.
    edge_tail : include the part of the edge average that lies below the grid,
           int_0^{v_min} H(g) dy, through the exact kernel of edge_tail_row() (exact for g
           supported on the grid).  DEFAULT True since 2026-11-22.
           With False (the historical matrices) J_op, forcing and Q_op equal the
           continuum operator MINUS (tail)*e^{-z}: Q_op's first row is identically 0 and
           the error is O(1) within a few units of z_min (for G = w/(2+w)^2, exact
           a = -0.17678: Q_op = 0.000 / -0.1115 / -0.1680 at z_min / z_min+1 / z_min+3).
           With True, Q_op G is accurate to ~1e-3 relative for BULK data.
           KNOWN REMAINING DEFECT (both settings): hilbert_support_z puts the p.v. of a
           constant over [1, Ymax] on the diagonal, i.e. row i assumes g = g_i on
           [1, Y_min], while the tail assumes g = 0 there.  On structure within O(1) of
           z_min Hg_op (hence Q_op) is 14-37% wrong, independent of N.  Edge-layer
           eigen/sup quantities (field of values, free-edge spectrum) are therefore
           grid-specific.  docs/findings_session_2026_11_22_qop_edge_tail.md section 6.
           Other effects of True: the stiff edge eigenvalue pair of M_GG grows ~1.5x
           (-141 +- 153i at z_min = -14).  Gates that add a tail themselves
           (Q_exact = Q_op + ...) must pass edge_tail=False.
    """
    z = np.linspace(zmin, zmax, N)
    dz = z[1] - z[0]
    Y = 1.0 + np.exp(z)
    ym1 = Y - 1.0                              # = e^z
    Ymax = Y[-1] + 0.5 * (Y[-1] - Y[-2])

    Dz = _dz_central(N, dz)
    Cz = _cumz(N, dz)
    Hz = hilbert_support_z(Y, Ymax)

    # g = (Y-1)^{-1/2} G ;  H(g) as a map on G
    Hg_op = Hz @ np.diag(ym1 ** -0.5)
    # J(Y) = int_1^Y H(g) dY' = int H(g) * (Y-1) dz'   (cumulative in z)
    J_op = Cz @ np.diag(ym1) @ Hg_op
    if edge_tail:                             # + int_0^{v_min} H(g) dy  (exact kernel)
        t_row = edge_tail_row(z)
        J_op = J_op + np.ones((N, 1)) @ t_row[None, :]
    forcing = np.diag(0.5 / ym1) @ J_op       # (1/2) J/(Y-1)

    # Bounded form of  delta_U/(Y-1) = J/(Y-1)  acting on G (= 2*forcing), recast as
    # the exponential-averaging convolution  int_{zmin}^{z} H(g) e^{z'-z} dz'  (weights
    # <= 1) so the matrix entries stay O(1) at the edge instead of being the discrete
    # small/small of diag(1/ym1)@J_op (cf. findings_transient_growth.md).  Used by the
    # nonlinear shape march for both the (verified) linear forcing and the nonlinear
    # self-advection term.  NOT beta-conjugated -> use at beta=0 only.
    EA = np.zeros((N, N))
    for i in range(1, N):
        w_ea = np.exp(z[:i + 1] - z[i])       # e^{z'-z} <= 1 for z' <= z
        tw = np.full(i + 1, dz); tw[0] *= 0.5; tw[i] *= 0.5   # trapezoid weights
        EA[i, :i + 1] = w_ea * tw
    Q_op = EA @ Hg_op                          # Q_op @ G ~= delta_U/(Y-1), bounded
    if edge_tail:
        Q_op = Q_op + np.outer(np.exp(-z), t_row)

    transport = -2.0 * Dz                     # -abar G_Y = -2 d_z (constant coeff)
    M_GG = transport + forcing
    M_hh = transport.copy()
    M_Gh = np.diag(ym1 ** -0.5) @ Dz          # +(Y-1)^{1/2} h_Y = e^{-z/2} d_z

    if beta != 0.0:                           # far-field weight similarity
        W = np.diag(np.exp(beta * z)); Wi = np.diag(np.exp(-beta * z))
        M_GG = W @ M_GG @ Wi; M_hh = W @ M_hh @ Wi; M_Gh = W @ M_Gh @ Wi

    if edge_bc:                               # G(z_min)=0, h(z_min)=0 (drop edge DOF)
        M_GG = M_GG[1:, 1:]; M_hh = M_hh[1:, 1:]; M_Gh = M_Gh[1:, 1:]
        Dz = Dz[1:, 1:]; Hg_op = Hg_op[1:, 1:]; J_op = J_op[1:, 1:]; Q_op = Q_op[1:, 1:]
        z = z[1:]; Y = Y[1:]

    return dict(M_GG=M_GG, M_hh=M_hh, M_Gh=M_Gh, z=z, Y=Y, dz=dz, beta=beta,
                Dz=Dz, Hg_op=Hg_op, J_op=J_op, Q_op=Q_op)


def spectrum(M):
    vals, vecs = np.linalg.eig(M)
    idx = np.argsort(-vals.real)
    return vals[idx], vecs[:, idx]


if __name__ == "__main__":
    print("edge-resolved (z) operator, inflow BC G(z_min)=0; far-field weight beta:")
    for beta in (0.0, -0.25, -0.5):
        B = build_blocks_z(N=1400, beta=beta)
        vGG, _ = spectrum(B["M_GG"])
        vhh, _ = spectrum(B["M_hh"])
        npos = int((vGG.real > 1e-6).sum())
        print(f"  beta={beta:+.2f}:  max Re eig(M_GG)={vGG.real.max():+.4f}  "
              f"(#Re>0={npos})   max Re eig(M_hh)={vhh.real.max():+.4f}   "
              f"edge/xi mode +{C_L-1:.0f}")
