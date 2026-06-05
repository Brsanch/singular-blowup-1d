#!/usr/bin/env python3
"""
Rescaled (self-similar) coordinates for the 1D Hou-Luo model, eq. (2.4) of
arXiv:2604.01868:

    Omega_tau + (U + c_l X) Omega_X = c_omega Omega + Theta_X
    Theta_tau + (U + c_l X) Theta_X = (c_l + 2 c_omega) Theta
    U_X = H(Omega),   U(0) = 0     (Hilbert transform on the real line R)

Steady states (Omega_tau = Theta_tau = 0) are the self-similar profiles.

This module provides:
  - line_hilbert : Hilbert transform on R for equispaced grids (validated against
    the Lorentzian, whose exact transform is x/(1+x^2)).
  - the EXPLICIT singular steady state of Theorem 2.3 and its analytic velocity,
  - steady_residual : residual of the steady eq. (2.4), to certify any profile.

Analytic facts (derived, see docs/findings_singular_steady_verified.md):
  Omega_bar(X) = 1{X>1} / sqrt(X-1),   Theta_bar(X) = (pi/2) 1{X>1},
  c_l = 2, c_omega = -1.
  U_X = H(Omega_bar) = -1/sqrt(1-X) for X<1,  0 for X>1.
  U(X) = 2 sqrt(1-X) - 2 for X<1,  -2 for X>1.   => U(1) = -2  (exactly the
  value the steady Theta-equation requires, since (c_l+2 c_omega)=0 forces
  (U+c_l X) Theta_bar_X = (U(1)+2)(pi/2) delta(X-1) = 0).
"""

import numpy as np


# ---------------------------------------------------------------------------
# Hilbert transform on the real line for an equispaced grid.
# H(f)(x) = (1/pi) p.v. \int f(y)/(x-y) dy.
# Equispaced "alternating" (odd-offset) quadrature is spectrally accurate for
# smooth, well-resolved, decaying f and is h-independent:
#     H(f)_j = (1/pi) * sum_{k != j} f_k / (j-k)   [trapezoidal p.v.]
# Using the odd-only rule  H(f)_j = (2/pi) sum_{j-k odd} f_k/(j-k)  is the
# standard accurate variant.  We build the (dense) kernel once; grids here are a
# few thousand points so O(N^2) is fine and exact about the p.v.
# ---------------------------------------------------------------------------
def hilbert_matrix(n):
    j = np.arange(n)
    M = np.zeros((n, n))
    d = j[:, None] - j[None, :]           # j - k
    odd = (d % 2 != 0)
    M[odd] = (2.0 / np.pi) / d[odd]
    return M


def line_hilbert(f, M=None):
    if M is None:
        M = hilbert_matrix(len(f))
    return M @ f


# ---------------------------------------------------------------------------
# The explicit singular steady state (Theorem 2.3) and its analytic velocity.
# eps regularizes the (X-1)^{-1/2} singularity for grid evaluation.
# ---------------------------------------------------------------------------
def singular_omega(X, eps=0.0):
    out = np.zeros_like(X)
    m = X > 1.0
    out[m] = 1.0 / np.sqrt(X[m] - 1.0 + eps)
    return out


def singular_theta(X):
    return np.where(X > 1.0, np.pi / 2.0, 0.0)


def singular_U_analytic(X):
    """U(X): 2 sqrt(1-X) - 2 for X<1, -2 for X>1."""
    return np.where(X < 1.0, 2.0 * np.sqrt(np.clip(1.0 - X, 0, None)) - 2.0, -2.0)


def singular_Ux_analytic(X):
    """U_X = H(Omega_bar): -1/sqrt(1-X) for X<1, 0 for X>1."""
    out = np.zeros_like(X)
    m = X < 1.0
    out[m] = -1.0 / np.sqrt(1.0 - X[m])
    return out


# ---------------------------------------------------------------------------
# Steady residual of eq. (2.4).  Pass analytic or numeric derivatives.
#   R_Omega = (U + c_l X) Omega_X - c_omega Omega - Theta_X
#   R_Theta = (U + c_l X) Theta_X - (c_l + 2 c_omega) Theta
# Evaluate away from the X=1 singularity (smooth region) to certify the profile.
# ---------------------------------------------------------------------------
def steady_residual(X, Omega, Omega_X, Theta, Theta_X, U, c_l, c_omega):
    drift = U + c_l * X
    R_Omega = drift * Omega_X - c_omega * Omega - Theta_X
    R_Theta = drift * Theta_X - (c_l + 2.0 * c_omega) * Theta
    return R_Omega, R_Theta
