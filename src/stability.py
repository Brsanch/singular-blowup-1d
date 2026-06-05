#!/usr/bin/env python3
"""
Linear-stability operator for the Hou-Luo singular self-similar profile
(arXiv:2604.01868, Thm 2.3 / Conjecture 2.4).

Linearize eq. (2.4) about the certified base state
    Omega_bar = 1{X>1}/sqrt(X-1),  Theta_bar = (pi/2) 1{X>1},
    U_bar(X)  = 2 sqrt(1-X) - 2 (X<1),  -2 (X>1)      [closed form, exact],
    c_l = 2,  c_omega = -1.
Frozen multipliers (mu = nu = 0; the scaling/normalization modes are added
separately).  With abar = U_bar + 2X:

    omega_tau = -abar omega_X - omega + theta_X - (H omega) Omega_bar_X
    theta_tau = -abar theta_X                   - (H omega) Theta_bar_X
    u_X = H(omega).   [u only enters through (H omega) here since u(0)=0 and the
                       coupling uses u' = H(omega'); we use H(omega) directly.]

Assembled as a 2N x 2N real matrix L on (omega, theta); eigenvalues with
Re(lambda) > 0 are unstable.  The base derivatives Omega_bar_X ~ (X-1)^{-3/2}
and Theta_bar_X ~ delta(X-1) are SINGULAR; we regularize with width-`w`/offset-
`delta` and study convergence (the genuine free-boundary subtlety lives here).

NOTE on u: linearizing (U + c_l X) with U=U_bar+u', the perturbation transport
term is (u' + mu X) Omega_bar_X; with mu=0 it is u' Omega_bar_X, u'=H(omega).
The base-velocity transport abar = U_bar + 2X is exact (closed form).
"""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.dirname(__file__))
from rescaled import hilbert_matrix  # noqa: E402


def abar(X):
    """Base transport speed U_bar + 2X (closed form)."""
    Ub = np.where(X < 1.0, 2.0 * np.sqrt(np.clip(1.0 - X, 0, None)) - 2.0, -2.0)
    return Ub + 2.0 * X


def _circ_deriv_central4(N, dx):
    """4th-order central first-derivative matrix (circulant)."""
    D = np.zeros((N, N))
    c = np.array([1.0, -8.0, 0.0, 8.0, -1.0]) / (12.0 * dx)
    offs = [-2, -1, 0, 1, 2]
    for a, o in zip(c, offs):
        D += a * np.eye(N, k=o)
        # wrap
        if o > 0:
            D += a * np.eye(N, k=o - N)
        elif o < 0:
            D += a * np.eye(N, k=o + N)
    return D


def _upwind3_matrix(N, dx, speed):
    """3rd-order upwind first-derivative matrix, branch per row by sign(speed)."""
    D = np.zeros((N, N))
    # a>=0 (left-biased): offsets {-2:1, -1:-6, 0:3, +1:2}/(6dx)
    pos = {-2: 1.0, -1: -6.0, 0: 3.0, 1: 2.0}
    # a<0 (right-biased): offsets {-1:-2, 0:-3, +1:6, +2:-1}/(6dx)
    neg = {-1: -2.0, 0: -3.0, 1: 6.0, 2: -1.0}
    for i in range(N):
        stencil = pos if speed[i] >= 0.0 else neg
        for o, val in stencil.items():
            j = (i + o) % N
            D[i, j] += val / (6.0 * dx)
    return D


def regularized_profile(X, delta=0.05, w=None):
    """Regularized base profile derivatives Omega_bar_X, Theta_bar_X.

    sig = smooth step at X=1 of width w; Omega = sig/sqrt(X-1+delta);
    Theta = (pi/2) sig.  Returns (Omega_X, Theta_X) on the grid via 4th-order
    central differences (consistent with the operator's source derivative).
    """
    if w is None:
        w = 2.0 * (X[1] - X[0])
    dx = X[1] - X[0]
    sig = 0.5 * (1.0 + np.tanh((X - 1.0) / w))
    denom = np.sqrt(np.clip(X - 1.0 + delta, delta / 50.0, None))
    Omega = sig / denom
    Theta = (np.pi / 2.0) * sig
    Dc = _circ_deriv_central4(len(X), dx)
    return Dc @ Omega, Dc @ Theta, Omega, Theta


def build_operator(N=1400, x0=-10.0, x1=60.0, delta=0.05, w=None, coupling=True):
    """Assemble the 2N x 2N linearized operator.  Returns (L, X)."""
    X = np.linspace(x0, x1, N)
    dx = X[1] - X[0]
    ab = abar(X)
    Dc = _circ_deriv_central4(N, dx)
    Dup = _upwind3_matrix(N, dx, ab)
    H = hilbert_matrix(N)
    I = np.eye(N)

    transport = -np.diag(ab) @ Dup           # -abar d_X
    A_ww = transport - I
    A_wt = Dc                                 # + theta_X
    A_tw = np.zeros((N, N))
    A_tt = transport.copy()

    if coupling:
        OmX, ThX, _, _ = regularized_profile(X, delta=delta, w=w)
        A_ww = A_ww - np.diag(OmX) @ H        # -(H omega) Omega_bar_X
        A_tw = A_tw - np.diag(ThX) @ H        # -(H omega) Theta_bar_X

    L = np.block([[A_ww, A_wt], [A_tw, A_tt]])
    return L, X


def spectrum(L, k=12):
    """Eigenvalues sorted by descending real part; return top-k with vectors."""
    vals, vecs = np.linalg.eig(L)
    idx = np.argsort(-vals.real)
    return vals[idx], vecs[:, idx]
