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


def build_blocks_z(N=1400, zmin=-14.0, zmax=7.0, beta=0.0, edge_bc=True):
    """Assemble the frozen operator on z in [zmin, zmax].

    beta : far-field weight exponent (w = e^{beta z}); beta<0 pushes the essential
           (transport) spectrum strictly left -> reveals the spectral gap.
    edge_bc : impose the inflow/regularity condition G(z_min)=0 (and h(z_min)=0).
           The edge z->-inf is an INFLOW boundary (characteristics move +z); the
           weighted space excludes the homogeneous branch G~(Y-1)^{-lambda/2} that
           blows up at the edge for Re lambda>0.  Without this BC the discretization
           re-admits that inadmissible branch as a spurious z_min-dependent edge
           mode (e.g. Re~+1.7).  Implemented by dropping the z_min DOF.
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
    forcing = np.diag(0.5 / ym1) @ J_op       # (1/2) J/(Y-1)

    transport = -2.0 * Dz                     # -abar G_Y = -2 d_z (constant coeff)
    M_GG = transport + forcing
    M_hh = transport.copy()
    M_Gh = np.diag(ym1 ** -0.5) @ Dz          # +(Y-1)^{1/2} h_Y = e^{-z/2} d_z

    if beta != 0.0:                           # far-field weight similarity
        W = np.diag(np.exp(beta * z)); Wi = np.diag(np.exp(-beta * z))
        M_GG = W @ M_GG @ Wi; M_hh = W @ M_hh @ Wi; M_Gh = W @ M_Gh @ Wi

    if edge_bc:                               # G(z_min)=0, h(z_min)=0 (drop edge DOF)
        M_GG = M_GG[1:, 1:]; M_hh = M_hh[1:, 1:]; M_Gh = M_Gh[1:, 1:]
        z = z[1:]; Y = Y[1:]

    return dict(M_GG=M_GG, M_hh=M_hh, M_Gh=M_Gh, z=z, Y=Y, dz=dz, beta=beta)


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
