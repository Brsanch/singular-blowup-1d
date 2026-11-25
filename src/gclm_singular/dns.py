#!/usr/bin/env python3
"""
Nonlinear dynamic-rescaling DNS solver for half-line gCLM (Case 2.2 of arXiv:2603.25104v4):

    Omega_tau + (c_l * X + a * U) * Omega_X = (c_omega + U_X) * Omega,
    U_X = H(Omega),    U(X, tau) = int_0^X H(Omega)(Y, tau) dY,

supported on X = 1 + v >= 1, with Case 2.2 normalization:
    c_l(tau) + a * U(1, tau) = 0,       H(Omega)(0, tau) = 1.

Discretization:
  * Exponentially stretched coordinate v(rho) = c_stretch * (e^rho - 1) on rho in [0, rhomax],
    resolving both the singular edge X = 1+ (v_1 ~ 1e-5) and the slow algebraic far-field tail
    out to v_max ~ 1e8.
  * 5th-order WENO (Jiang-Shu) for advection respecting the local sign of A(1+v) = c_l*v + a*J(v),
    coupled with TVD 3-stage SSP-RK3 in time.
  * Exact log-kernel quadrature for J(v) = U(1+v) - U(1) = (1/pi) int Om(1+v') ln|(v-v')/v'| dv'
    (avoids evaluating the divergent H(Omega)(1) at the singular edge).
"""

import numpy as np


def weno5_deriv(f: np.ndarray, speed: np.ndarray, h: float) -> np.ndarray:
    """5th-order WENO first derivative f_rho respecting the local sign of speed."""
    N = len(f)
    fp = np.empty(N + 6)
    fp[3:-3] = f
    fp[:3] = 0.0
    fp[-3:] = f[-1]
    df = np.diff(fp) / h
    eps = 1e-12

    def wm(v1, v2, v3, v4, v5):
        p1 = (1.0 / 3.0) * v1 - (7.0 / 6.0) * v2 + (11.0 / 6.0) * v3
        p2 = -(1.0 / 6.0) * v2 + (5.0 / 6.0) * v3 + (1.0 / 3.0) * v4
        p3 = (1.0 / 3.0) * v3 + (5.0 / 6.0) * v4 - (1.0 / 6.0) * v5
        s1 = (13.0 / 12.0) * (v1 - 2.0 * v2 + v3) ** 2 + 0.25 * (v1 - 4.0 * v2 + 3.0 * v3) ** 2
        s2 = (13.0 / 12.0) * (v2 - 2.0 * v3 + v4) ** 2 + 0.25 * (v2 - v4) ** 2
        s3 = (13.0 / 12.0) * (v3 - 2.0 * v4 + v5) ** 2 + 0.25 * (3.0 * v3 - 4.0 * v4 + v5) ** 2
        a1 = 0.1 / (eps + s1) ** 2
        a2 = 0.6 / (eps + s2) ** 2
        a3 = 0.3 / (eps + s3) ** 2
        return (a1 * p1 + a2 * p2 + a3 * p3) / (a1 + a2 + a3)

    dm = wm(df[0:N], df[1 : N + 1], df[2 : N + 2], df[3 : N + 3], df[4 : N + 4])
    dp = wm(df[5 : N + 5], df[4 : N + 4], df[3 : N + 3], df[2 : N + 2], df[1 : N + 1])
    return np.where(speed >= 0.0, dm, dp)


def build_exact_log_J(vi: np.ndarray, wv: np.ndarray) -> np.ndarray:
    """Exact log-kernel matrix for J(v_i) = (1/pi) int Om(1+v') ln|(v_i - v')/v'| dv'."""
    n = len(vi)
    diff = np.abs(vi[:, None] - vi[None, :])
    with np.errstate(divide="ignore", invalid="ignore"):
        L = np.log(diff / vi[None, :])
    # Exact diagonal cell integral of ln|v_i - v'| - ln(v_i) over width h_i = wv[i]
    L[np.diag_indices(n)] = np.log(0.5 * wv / vi) - 1.0
    return (L * wv[None, :]) / np.pi


def run_case22_dns(
    a: float,
    init_type: str = "exp_inv",
    N: int = 450,
    c_stretch: float = 1e-4,
    rhomax: float = 29.0,
    T_end: float = 15.0,
    cfl: float = 0.25,
) -> dict:
    """Evolve Case 2.2 nonlinear dynamic-rescaling gCLM on v(rho) = c_stretch * (e^rho - 1).

    init_type:
        'exp_inv'  : arXiv:2603.25104 Fig. 2(d) Case 2.2 datum -X^{-3} exp(-1/(X-1))
        'poly_exp' : smooth degenerate polynomial-exponential -(X-1)^4 exp(-(X-1))
        'branch2'  : regularized Branch II datum -(X - 1 + 1e-3)^{-1/2}
    """
    a = float(a)
    rho = np.linspace(0.0, rhomax, N)
    dr = rho[1] - rho[0]
    v = c_stretch * np.expm1(rho)
    dv = c_stretch * np.exp(rho)
    X = 1.0 + v
    v_right = v[-1] + 0.5 * (v[-1] - v[-2])

    vi = v[1:]
    dvi = dv[1:]
    Xi = X[1:]
    ni = N - 1
    w_r = np.full(ni, dr)
    w_r[-1] *= 0.5
    wv = dvi * w_r

    # Subtraction-form Hilbert transform on [0, v_right]
    d = vi[:, None] - vi[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        K = wv[None, :] / d
    K[np.diag_indices(ni)] = 0.0
    Hv = K / np.pi
    Hv[np.diag_indices(ni)] = (np.log(vi / (v_right - vi)) - K.sum(axis=1)) / np.pi

    J_op = build_exact_log_J(vi, wv)
    u1_row = -(1.0 / np.pi) * np.log1p(1.0 / vi) * wv
    l0_row = -(1.0 / np.pi) * (wv / (1.0 + vi))

    if init_type == "exp_inv":
        Om = -(Xi ** -3) * np.exp(-1.0 / np.maximum(vi, 1e-300))
    elif init_type == "poly_exp":
        Om = -(vi ** 4) * np.exp(-vi)
    elif init_type == "branch2":
        Om = -((vi + 1e-3) ** -0.5)
    else:
        raise ValueError(f"Unknown init_type={init_type!r}")
    Om /= float(l0_row @ Om)

    def rhs(O):
        Ux = Hv @ O
        U1 = float(u1_row @ O)
        J = J_op @ O
        c_l = -a * U1
        A = c_l * vi + a * J
        spd = A / dvi
        r0 = -spd * weno5_deriv(O, spd, dr) + Ux * O
        c_w = -float(l0_row @ r0) / float(l0_row @ O)
        return r0 + c_w * O, c_l, c_w, A, Ux

    tau = 0.0
    steps = 0
    while tau < T_end:
        k1, c_l, c_w, A, Ux = rhs(Om)
        dt = min(
            cfl * dr / (float(np.max(np.abs(A / dvi))) + 1e-12),
            0.12 / (float(np.max(np.abs(c_w + Ux))) + 1e-12),
            0.02,
            T_end - tau,
        )
        O1 = Om + dt * k1
        k2, *_ = rhs(O1)
        O2 = 0.75 * Om + 0.25 * (O1 + dt * k2)
        k3, *_ = rhs(O2)
        Om = (1.0 / 3.0) * Om + (2.0 / 3.0) * (O2 + dt * k3)
        Om /= float(l0_row @ Om)
        tau += dt
        steps += 1
        if not np.all(np.isfinite(Om)) or steps > 60000:
            break

    r_fin, c_l, c_w, A, Ux = rhs(Om)
    gamma = -c_l / c_w

    # Power-law slope in bulk v in [0.2, 30]
    m_fit = (vi > 0.2) & (vi < 30.0) & (Om < -1e-10)
    beta_fit = (
        float(-np.polyfit(np.log(vi[m_fit]), np.log(-Om[m_fit]), 1)[0])
        if m_fit.sum() > 5
        else np.nan
    )

    # Interior stagnation point X_* where A changes sign from <0 to >=0
    sc = np.where((A[:-1] < 0.0) & (A[1:] >= 0.0))[0]
    x_star = float(1.0 + vi[sc[-1]]) if len(sc) > 0 else np.nan

    m_bulk = (vi > 0.05) & (vi < 50.0)
    res_bulk = float(np.max(np.abs(r_fin[m_bulk])))

    return dict(
        a=a, init_type=init_type, tau=tau, steps=steps,
        v=vi, X=Xi, Om=Om, A=A, Ux=Ux,
        c_l=c_l, c_w=c_w, gamma=gamma, beta_fit=beta_fit,
        x_star=x_star, res_bulk=res_bulk,
    )
