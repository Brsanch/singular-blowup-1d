#!/usr/bin/env python3
"""
Exact singular self-similar profiles for the generalized Constantin-Lax-Majda
(gCLM) model (Huang-Tong-Wang, arXiv:2603.25104v4):

    omega_t + a * u * omega_x = u_x * omega,    u_x = H(omega),    u(0,t) = 0,

in dynamic-rescaling coordinates:

    Omega_tau + (c_l * X + a * U) * Omega_X = (c_omega + U_X) * Omega,
    U_X = H(Omega),    U(X,tau) = int_0^X H(Omega)(Y,tau) dY.

Hilbert transform convention (identical to arXiv:2603.25104 and src/gclm.py):
    H(f)(X) = (1/pi) p.v. int_{-inf}^{inf} f(Y) / (X - Y) dY.

Theorem 3 (half-line explicit singular family, a < 0):
    mu = 1 / (1 - a) in (0, 1),
    Omega_bar_a(X) = -sin(pi * mu) * (X - 1)^{-mu} * 1_{X > 1},
    c_l_bar = 1 - a,    c_omega_bar = -1,    gamma_bar = -c_l_bar / c_omega_bar = 1 - a.
"""

import numpy as np
from scipy.integrate import quad


def mu_of_a(a: float) -> float:
    """Singularity exponent mu(a) = 1/(1-a) in (0, 1) for a < 0."""
    a = float(a)
    if a >= 0.0:
        raise ValueError(f"Explicit half-line singular profile requires a < 0; got a={a}")
    return 1.0 / (1.0 - a)


def scaling_exponents(a: float):
    """Return (c_l, c_omega, gamma) = (1 - a, -1, 1 - a) for a < 0."""
    a = float(a)
    if a >= 0.0:
        raise ValueError(f"Expected a < 0; got a={a}")
    c_l = 1.0 - a
    c_omega = -1.0
    gamma = -c_l / c_omega
    return c_l, c_omega, gamma


def omega_bar(X, a: float):
    """Exact singular profile Omega_bar_a(X) = -sin(pi*mu) * (X-1)^{-mu} * 1_{X>1}."""
    X = np.asarray(X, dtype=float)
    mu = mu_of_a(a)
    out = np.zeros_like(X, dtype=float)
    mask = X > 1.0
    out[mask] = -np.sin(np.pi * mu) * (X[mask] - 1.0) ** (-mu)
    return out


def omega_bar_X(X, a: float):
    """Exact derivative (Omega_bar_a)_X(X) = mu*sin(pi*mu)*(X-1)^{-1-mu}*1_{X>1}."""
    X = np.asarray(X, dtype=float)
    mu = mu_of_a(a)
    out = np.zeros_like(X, dtype=float)
    mask = X > 1.0
    out[mask] = mu * np.sin(np.pi * mu) * (X[mask] - 1.0) ** (-1.0 - mu)
    return out


def Ux_bar(X, a: float):
    """Exact Hilbert transform H(Omega_bar_a)(X) = (1_{X<1} + cos(pi*mu)*1_{X>1}) / |1-X|^mu."""
    X = np.asarray(X, dtype=float)
    mu = mu_of_a(a)
    out = np.empty_like(X, dtype=float)
    left = X < 1.0
    right = X > 1.0
    out[left] = (1.0 - X[left]) ** (-mu)
    out[right] = np.cos(np.pi * mu) * (X[right] - 1.0) ** (-mu)
    out[~(left | right)] = np.nan
    return out


def U_bar(X, a: float):
    """Exact velocity U_bar_a(X) = int_0^X H(Omega_bar_a)(Y) dY, U_bar_a(0) = 0."""
    X = np.asarray(X, dtype=float)
    mu = mu_of_a(a)
    om_mu = 1.0 - mu
    out = np.empty_like(X, dtype=float)
    left = X <= 1.0
    right = ~left
    out[left] = (1.0 - (1.0 - X[left]) ** om_mu) / om_mu
    out[right] = (1.0 + np.cos(np.pi * mu) * (X[right] - 1.0) ** om_mu) / om_mu
    return out


def advection_bar(X, a: float):
    """Background advection speed A_bar_a(X) = c_l_bar * X + a * U_bar_a(X).

    On X = 1 + v > 1, equals (1 - a) * (v - cos(pi*mu) * v^{1 - mu}).
    """
    c_l, _, _ = scaling_exponents(a)
    X = np.asarray(X, dtype=float)
    return c_l * X + a * U_bar(X, a)


def steady_residual(X, a: float):
    """Pointwise steady residual R(X) = -(c_l*X + a*U)*Omega_X + (c_omega + U_X)*Omega."""
    c_l, c_omega, _ = scaling_exponents(a)
    A = advection_bar(X, a)
    Om = omega_bar(X, a)
    OmX = omega_bar_X(X, a)
    Ux = Ux_bar(X, a)
    return -A * OmX + (c_omega + Ux) * Om


def interior_stagnation_point(a: float):
    """Return the interior stagnation point X_* > 1 of A_bar_a(X) for a < -1.

    For -1 < a < 0: cos(pi*mu) < 0, A_bar_a > 0 on all (1, inf) -> returns None.
    For a == -1:    cos(pi*mu) = 0, A_bar_{-1}(X) = 2*(X-1) -> returns 1.0 (edge only).
    For a < -1:     cos(pi*mu) > 0, A_bar_a vanishes at X_* = 1 + (cos(pi*mu))^{1/mu} > 1.
    """
    mu = mu_of_a(a)
    c = np.cos(np.pi * mu)
    if np.isclose(c, 0.0, atol=1e-15):
        return 1.0
    if c < 0.0:
        return None
    return 1.0 + c ** (1.0 / mu)


def hilbert_omega_bar_quad(X_eval, a: float):
    """Independent high-precision adaptive-quadrature evaluation of H(Omega_bar_a)(X_eval).

    Uses the substitution v = t^p with p = 1/(1-mu) on the singular edge [0, v_split]
    to completely eliminate the integrable v^{-mu} singularity, and u = 1/v on [v_far, inf)
    to capture the algebraic tail without truncation error.
    For X_eval > 1 (principal value at v0 = X_eval - 1), uses symmetric subtraction around v0.
    """
    mu = mu_of_a(a)
    pref = -np.sin(np.pi * mu) / np.pi
    X_arr = np.atleast_1d(np.asarray(X_eval, dtype=float))
    res = np.empty_like(X_arr)
    p = 1.0 / (1.0 - mu)  # v = t^p -> v^{-mu} dv = p * dt (completely smooth at t=0!)

    for idx, x in enumerate(X_arr):
        if x < 1.0:
            # Integral = pref * int_0^inf v^{-mu} / (x - 1 - v) dv
            # Split at v = 1: on [0, 1] set v = t^p (t in [0,1]); on [1, inf) set v = 1/u (u in (0,1])
            d = 1.0 - x  # > 0, denominator is -(d + v)
            I1, _ = quad(lambda t: -p / (d + t ** p), 0.0, 1.0, epsabs=1e-13, epsrel=1e-13)
            # v = 1/u -> dv = -du/u^2, v^{-mu}/(x-1-v) dv = - u^{mu-2} / (d + 1/u) du = - u^{mu-1}/(d*u + 1) du
            # Substitute u = s^q with q = 1/mu so u^{mu-1} du = q ds!
            q = 1.0 / mu
            I2, _ = quad(lambda s: -q / (d * (s ** q) + 1.0), 0.0, 1.0, epsabs=1e-13, epsrel=1e-13)
            res[idx] = pref * (I1 + I2)
        elif x > 1.0:
            v0 = x - 1.0  # > 0
            # p.v. int_0^inf v^{-mu} / (v0 - v) dv
            # Split into:
            #   (a) [0, v0/2]: v = (v0/2) * t^p
            #   (b) [v0/2, 3*v0/2]: symmetric p.v. via h in [0, v0/2]:
            #       int_0^{v0/2} ((v0-h)^{-mu} - (v0+h)^{-mu}) / h dh  (smooth at h=0!)
            #   (c) [3*v0/2, inf): v = 1/u, u in [0, 2/(3*v0)], then u = (2/(3*v0)) * s^q
            a1 = 0.5 * v0
            I1, _ = quad(
                lambda t: p * (a1 ** (1.0 - mu)) / (v0 - a1 * (t ** p)),
                0.0, 1.0, epsabs=1e-13, epsrel=1e-13,
            )
            I2, _ = quad(
                lambda h: ((v0 - h) ** (-mu) - (v0 + h) ** (-mu)) / h,
                0.0, a1, epsabs=1e-13, epsrel=1e-13,
            )
            u_max = 1.0 / (1.5 * v0)
            q = 1.0 / mu
            # int_{3v0/2}^inf v^{-mu}/(v0 - v) dv = int_0^{u_max} u^{mu-1}/(v0*u - 1) du
            I3, _ = quad(
                lambda s: q * (u_max ** mu) / (v0 * u_max * (s ** q) - 1.0),
                0.0, 1.0, epsabs=1e-13, epsrel=1e-13,
            )
            res[idx] = pref * (I1 + I2 + I3)
        else:
            res[idx] = np.nan

    return float(res[0]) if np.ndim(X_eval) == 0 else res


def mellin_symbol_c(s, a: float):
    """Exact Mellin shift symbol c_a(s) of M_a^{bulk} for 0 < Re(s) < 1:

        M[M_a^{bulk} G](s) = (s / mu) * G_hat(s) + c_a(s) * G_hat(s - mu),

        c_a(s) = (mu - s) * [ cos(pi*mu)/mu + sin(pi*mu) * cot(pi*s) / (1 - s) ].
    """
    mu = mu_of_a(a)
    s = np.asarray(s, dtype=complex)
    cot_pi_s = np.cos(np.pi * s) / np.sin(np.pi * s)
    return (mu - s) * (np.cos(np.pi * mu) / mu + np.sin(np.pi * mu) * cot_pi_s / (1.0 - s))


def conjugation_N_am1(s):
    """Exact elementary Mellin conjugation multiplier at a = -1 (mu = 1/2):

        N_{-1}(s) = sin(pi*s) / (2*s - 1),

    satisfying N_{-1}(s) * c_{-1}(s) = -N_{-1}(s - 1/2).
    """
    s = np.asarray(s, dtype=complex)
    return np.sin(np.pi * s) / (2.0 * s - 1.0)
