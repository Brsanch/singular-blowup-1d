#!/usr/bin/env python3
"""
Desingularized free-boundary stability operator for the Hou-Luo singular profile
(arXiv:2604.01868, Conjecture 2.4).  Built from the corrected edge-condition
derivation in docs/findings_edge_condition_derivation.md (verified by
validate/probe_edge_condition.py).

KEY STRUCTURE (why this is well-posed where the naive operator was not):

1. MOVING (edge-fixed) frame.  The fixed-X-frame edge shift -xi*Omega_bar' has a
   SPURIOUS infinite edge velocity (Taylor breakdown at the vertical-tangent edge)
   -> that was the delta-divergent edge artifact.  In the moving frame the
   perturbation velocity is finite.

2. xi DROPS OUT of the shape equation.  Substituting the edge ODE
   xi_dot = (c_l-1) xi + V_g(1) + mu  into the Omega-equation turns the singular
   edge term into  -[ int_1^Y H(g) dY' ] Omega_bar_Y , a pure functional of g.
   The operator is BLOCK-TRIANGULAR:  the (g,h) block is self-contained, and the
   full spectrum is  eig(M) U {c_l-1}.  The pure-edge (xi) symmetry mode at
   lambda = c_l-1 = 1 is appended analytically (NOT discretized) -> the edge
   singularity is never represented on the grid.

3. WEIGHTED variable  G = (Y-1)^{1/2} g  (so g = (Y-1)^{-1/2} G).  Because
   c_omega = -1 exactly cancels the weight's contribution, the G-equation has NO
   zeroth-order term and a BOUNDED forcing:
       G_tau = -abar G_Y + (Y-1)^{1/2} h_Y + (1/2) J(Y)/(Y-1),   J=int_1^Y H(g),
   with (1/2) J/(Y-1) -> (1/2) H(g)(1) at the edge (finite).  h (the Theta-shape)
   is regular at the edge, unweighted:  h_tau = -abar h_Y  (frozen mu=nu=0).

Frozen multipliers mu = nu = 0 here (the scaling/normalization modes mu, nu are a
separate, later addition).  Far field: NON-periodic outflow derivative (abar>0 on
the support, outgoing at Y=L) -> no circulant wrap, so NO spurious box modes.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from rescaled import hilbert_matrix  # noqa: E402

C_L = 2.0
C_OMEGA = -1.0


def _backward_upwind3_matrix(N, dx):
    """3rd-order backward-biased d_Y on a NON-periodic grid (abar>0 everywhere on
    the support: outgoing at the right end -> one-sided there; abar->0 at the left
    edge so its small error is harmless).  No circulant wrap."""
    D = np.zeros((N, N))
    # interior i>=2 and i<=N-2 : (2 f_{i+1} + 3 f_i - 6 f_{i-1} + f_{i-2})/(6 dx)
    for i in range(2, N - 1):
        D[i, i + 1] += 2.0 / (6 * dx)
        D[i, i] += 3.0 / (6 * dx)
        D[i, i - 1] += -6.0 / (6 * dx)
        D[i, i - 2] += 1.0 / (6 * dx)
    # i = 0 (edge side): 2nd-order forward (abar ~ 0 here, weight tiny)
    D[0, 0] += -3.0 / (2 * dx)
    D[0, 1] += 4.0 / (2 * dx)
    D[0, 2] += -1.0 / (2 * dx)
    # i = 1 : 2nd-order central
    D[1, 2] += 1.0 / (2 * dx)
    D[1, 0] += -1.0 / (2 * dx)
    # i = N-1 (outflow right end): 2nd-order backward one-sided
    D[N - 1, N - 1] += 3.0 / (2 * dx)
    D[N - 1, N - 2] += -4.0 / (2 * dx)
    D[N - 1, N - 3] += 1.0 / (2 * dx)
    return D


def _central2_matrix(N, dx):
    """2nd-order central first-derivative, non-periodic (one-sided at ends)."""
    D = np.zeros((N, N))
    for i in range(1, N - 1):
        D[i, i + 1] += 1.0 / (2 * dx)
        D[i, i - 1] += -1.0 / (2 * dx)
    D[0, 0] += -3.0 / (2 * dx); D[0, 1] += 4.0 / (2 * dx); D[0, 2] += -1.0 / (2 * dx)
    D[N - 1, N - 1] += 3.0 / (2 * dx); D[N - 1, N - 2] += -4.0 / (2 * dx)
    D[N - 1, N - 3] += 1.0 / (2 * dx)
    return D


def _cumtrapz_matrix(N, dx, Y):
    """Matrix C such that (C @ f)_i ~ int_1^{Y_i} f dY', with the first sub-cell
    [1, Y_0] approximated by the rectangle f_0*(Y_0-1) (Y_0 = 1 + h0)."""
    C = np.zeros((N, N))
    # J_0 = f_0 * (Y_0 - 1)   (rectangle from the edge Y=1 to the first node)
    C[0, 0] = (Y[0] - 1.0)
    # J_i = J_{i-1} + 0.5 (f_i + f_{i-1}) dx
    for i in range(1, N):
        C[i] = C[i - 1]
        C[i, i] += 0.5 * dx
        C[i, i - 1] += 0.5 * dx
    return C


def build_blocks(N=1400, h0=2e-3, L=120.0):
    """Assemble the frozen (mu=nu=0) shape operator on the support Y in [1+h0, L].

    Returns dict with M_GG (G-block, weighted), M_hh (h-block), M_Gh (h->G source),
    and the grid Y.  Eigenvalues of the full frozen operator are
        eig(M_GG) U eig(M_hh) U {C_L - 1}.
    The G->h coupling is zero (block-triangular), so eigenVALUES need only the
    diagonal blocks + the appended edge mode; M_Gh is returned for eigenvectors.
    """
    Y = np.linspace(1.0 + h0, L, N)
    dx = Y[1] - Y[0]
    ym1 = Y - 1.0                       # (Y-1) > 0 on the support
    abar = 2.0 * ym1                    # closed form on Y>1: U_bar + 2Y = 2(Y-1)

    Dup = _backward_upwind3_matrix(N, dx)
    Dc = _central2_matrix(N, dx)
    Hful = hilbert_matrix(N)            # line-Hilbert on this (support-only) grid
    Cint = _cumtrapz_matrix(N, dx, Y)
    I = np.eye(N)

    # --- M_GG : G_tau = -abar G_Y + (1/2) J(Y)/(Y-1),  J = int_1^Y H(g),
    #            g = (Y-1)^{-1/2} G.  All operators in terms of G.
    W_inv = np.diag(ym1 ** -0.5)        # g = W_inv @ G   (g = (Y-1)^{-1/2} G)
    Hg_op = Hful @ W_inv               # H(g) as a linear map on G
    J_op = Cint @ Hg_op                # J(Y) = int_1^Y H(g)  as a map on G
    forcing_GG = np.diag(0.5 / ym1) @ J_op     # (1/2) J/(Y-1)  (bounded at edge)
    transport_GG = -np.diag(abar) @ Dup
    M_GG = transport_GG + forcing_GG

    # --- M_hh : h_tau = -abar h_Y   (Theta-shape, unweighted, regular at edge)
    M_hh = -np.diag(abar) @ Dup

    # --- M_Gh : + (Y-1)^{1/2} h_Y   (h sources G; not needed for eigenVALUES)
    M_Gh = np.diag(ym1 ** 0.5) @ Dc

    return dict(M_GG=M_GG, M_hh=M_hh, M_Gh=M_Gh, Y=Y, dx=dx, abar=abar,
                W_inv=W_inv, Hg_op=Hg_op, J_op=J_op)


def spectrum(M, k=None):
    vals, vecs = np.linalg.eig(M)
    idx = np.argsort(-vals.real)
    if k is not None:
        idx = idx[:k]
    return vals[idx], vecs[:, idx]


def full_frozen_spectrum(N=1400, h0=2e-3, L=120.0):
    """Frozen (mu=nu=0) spectrum of the desingularized operator:
       eig(M_GG) U eig(M_hh) U {C_L - 1 (edge/xi symmetry mode)}.
    Returns (all_vals_sorted, blocks)."""
    B = build_blocks(N=N, h0=h0, L=L)
    vGG, _ = spectrum(B["M_GG"])
    vhh, _ = spectrum(B["M_hh"])
    edge_mode = np.array([C_L - 1.0 + 0j])      # the appended pure-xi mode
    allv = np.concatenate([vGG, vhh, edge_mode])
    allv = allv[np.argsort(-allv.real)]
    return allv, B


if __name__ == "__main__":
    allv, B = full_frozen_spectrum()
    print("desingularized frozen spectrum (top 8 by Re):")
    for v in allv[:8]:
        print(f"   {v.real:+.5f} {v.imag:+.5f}i")
    vGG, _ = spectrum(B["M_GG"])
    print(f"\n  max Re eig(M_GG)  = {vGG.real.max():+.5f}")
    print(f"  appended edge mode = {C_L-1:+.5f}  (the symmetry mode to project out)")
