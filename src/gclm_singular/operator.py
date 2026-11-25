#!/usr/bin/env python3
"""
Desingularized linearized stability operator M_a for the explicit singular
self-similar profiles of the gCLM model (arXiv:2603.25104v4, Theorem 3, a < 0).

Support-interior perturbation on X = 1 + v > 1:
    Omega(1+v, tau) = Omega_bar_a(1+v) * (1 + G(1+v, tau)),
    g(1+v, tau)     = Omega_bar_a(1+v) * G(1+v, tau) = -sin(pi*mu) * v^{-mu} * G(1+v, tau),
    mu = 1 / (1 - a) in (0, 1),   a * mu = -(1 - mu).

Under edge pinning c_l + a*U(1) = 0 (delta_c_l = -a * delta_U(1)):
    G_tau = M_a^{bulk} G + kappa[G],
    M_a^{bulk} G(1+v) = -A_bar_a(1+v) G_v + H_+[g](1+v) + a*mu * Q_+[g](1+v),
    A_bar_a(1+v) = (1 - a) * (v - cos(pi*mu) * v^{1 - mu}),
    Q_+[g](1+v)  = (1 / v) * int_0^v H_+[g](1+y) dy,
    kappa[G]     = delta_c_omega - a*mu * delta_U(1).

Gauges:
    gauge='a' (gCLM analogue of Hou-Luo Theta-balance):
        c_omega = -mu * c_l  =>  kappa == 0 identically (even nonlinearly).
        M_a = M_a^{bulk}, with amplitude scaling mode G == 1 in ker(M_a).
    gauge='b' (arXiv:2603.25104 Case 2.2 normalization):
        c_l + a*U(1) = 0,  H(Omega)(0) == 1  =>  L_0[g] = 0,
        kappa[G] = -L_0[Omega_bar_a * M_a^{bulk} G]  (rank-1 projection along 1).

Coordinates:
    1. z-grid:    z = ln(X - 1) in [zmin, zmax]  (works for ALL a < 0).
                  Transport: -(1-a)*(1 - cos(pi*mu)*e^{-mu z}) * d_z.
                  At a = -1: pure constant-speed -2*d_z (logarithmic inflow at z -> -inf).
                  For -1 < a < 0: strictly positive speed (inflow at z_min).
                  For a < -1: changes sign at z_* = (1/mu)*ln(cos(pi*mu)) (outflow at BOTH ends!).
    2. zeta-grid: zeta = ln(1 + v^mu / |cos(pi*mu)|) in [0, zetamax]  (for -1 < a < 0).
                  Exact characteristic coordinate: transport is literally -d_zeta (speed 1),
                  with exact physical edge X = 1 at zeta = 0 (inflow BC G(0) = 0).
"""

import numpy as np
from scipy.integrate import quad

from .profiles import mu_of_a, scaling_exponents, omega_bar, omega_bar_X, U_bar, Ux_bar, advection_bar


def _d1_central(N: int, h: float) -> np.ndarray:
    """2nd-order central first derivative on a uniform grid, 2nd-order one-sided at ends."""
    D = np.zeros((N, N))
    for i in range(1, N - 1):
        D[i, i + 1] = 0.5 / h
        D[i, i - 1] = -0.5 / h
    D[0, 0] = -1.5 / h
    D[0, 1] = 2.0 / h
    D[0, 2] = -0.5 / h
    D[-1, -1] = 1.5 / h
    D[-1, -2] = -2.0 / h
    D[-1, -3] = 0.5 / h
    return D


def _d1_upwind2(N: int, h: float, speed: np.ndarray) -> np.ndarray:
    """2nd-order upwind first derivative matching the sign of advection speed c(x) in -c(x) d_x.

    For speed[i] >= 0 (characteristics moving +x), uses backward stencil (3f_i - 4f_{i-1} + f_{i-2})/(2h).
    For speed[i] < 0  (characteristics moving -x), uses forward stencil  (-3f_i + 4f_{i+1} - f_{i+2})/(2h).
    Falls back to 1st-order upwind at the 1-cell inflow margin.
    """
    D = np.zeros((N, N))
    for i in range(N):
        if speed[i] >= 0.0:
            if i >= 2:
                D[i, i] = 1.5 / h
                D[i, i - 1] = -2.0 / h
                D[i, i - 2] = 0.5 / h
            elif i == 1:
                D[i, 1] = 1.0 / h
                D[i, 0] = -1.0 / h
            else:
                D[0, 0] = -1.5 / h
                D[0, 1] = 2.0 / h
                D[0, 2] = -0.5 / h
        else:
            if i <= N - 3:
                D[i, i] = -1.5 / h
                D[i, i + 1] = 2.0 / h
                D[i, i + 2] = -0.5 / h
            elif i == N - 2:
                D[i, N - 2] = -1.0 / h
                D[i, N - 1] = 1.0 / h
            else:
                D[-1, -1] = 1.5 / h
                D[-1, -2] = -2.0 / h
                D[-1, -3] = 0.5 / h
    return D


def hilbert_support_nonuniform(v: np.ndarray, dv_dxi: np.ndarray, dxi: float,
                               v_left: float, v_right: float,
                               grid_support: bool = True) -> np.ndarray:
    """Subtraction-form line Hilbert transform for g supported on v in [v_left, v_right].

    Evaluates H_+[g](1 + v_i) = (1/pi) p.v. int g(1+v') / (v_i - v') dv'
    using exact coordinate Jacobian dv = dv_dxi * dxi and uniform trapezoid weights.
    Exact on constants over [v_left, v_right] (or [v[0], v_right] when grid_support=True).
    """
    n = len(v)
    w_xi = np.full(n, dxi)
    w_xi[0] *= 0.5
    w_xi[-1] *= 0.5
    w_v = dv_dxi * w_xi

    d = v[:, None] - v[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        K = w_v[None, :] / d
    K[np.diag_indices(n)] = 0.0

    # Exact p.v. of constant 1 on [va, v_right]: ln((v_i - va) / (v_right - v_i))
    pv = np.empty(n)
    if grid_support:
        # On row 0, v[0] - v[0] = 0 is log-divergent; use v_left (0.0) if v[0] > v_left, else 0.0
        pv[0] = np.log((v[0] - v_left) / (v_right - v[0])) if v[0] > v_left else 0.0
        pv[1:] = np.log((v[1:] - v[0]) / (v_right - v[1:]))
    else:
        mask = v > v_left
        pv[mask] = np.log((v[mask] - v_left) / (v_right - v[mask]))
        pv[~mask] = 0.0

    H = K / np.pi
    H[np.diag_indices(n)] = (pv - K.sum(axis=1)) / np.pi
    return H


def edge_tail_row_z(z: np.ndarray, mu: float) -> np.ndarray:
    """Row vector t such that t @ G = int_0^{v_0} H_+[g](1 + y) dy (v_0 = e^{z_0}),
    EXACT for g(1+v) = -sin(pi*mu) * v^{-mu} * G(1+v) supported on v >= v_0:

        t @ G = +(sin(pi*mu) / pi) * int_{z_0}^{z_max} e^{(1-mu)*z'} * (-ln(1 - e^{z_0 - z'})) * G(z') dz'.

    Piecewise-linear hat-function weights (12-point Gauss-Legendre on smooth cells,
    adaptive quad on the singular first cell [z_0, z_1]).
    """
    z = np.asarray(z, dtype=float)
    n = len(z)
    z0 = z[0]
    om_mu = 1.0 - mu

    def kernel(zp):
        return np.exp(om_mu * zp) * (-np.log1p(-np.exp(z0 - zp)))

    xg, wg = np.polynomial.legendre.leggauss(12)
    W = np.zeros(n)

    # Smooth cells j = 1 .. n-2
    a_c, b_c = z[1:-1], z[2:]
    half = 0.5 * (b_c - a_c)
    zq = a_c[:, None] + half[:, None] * (xg[None, :] + 1.0)
    kq = kernel(zq) * wg[None, :] * half[:, None]
    lam = (zq - a_c[:, None]) / (b_c - a_c)[:, None]
    W[1:-1] += (kq * (1.0 - lam)).sum(axis=1)
    W[2:] += (kq * lam).sum(axis=1)

    # Singular first cell [z_0, z_1] (integrable log singularity at z_0)
    h = z[1] - z[0]
    W[0] += quad(lambda zp: kernel(zp) * (z[1] - zp) / h, z0, z[1], limit=200)[0]
    W[1] += quad(lambda zp: kernel(zp) * (zp - z0) / h, z0, z[1], limit=200)[0]

    return (np.sin(np.pi * mu) / np.pi) * W


def delta_U1_row(v: np.ndarray, dv_dxi: np.ndarray, dxi: float, mu: float) -> np.ndarray:
    """Row vector u1 such that u1 @ G = delta_U(1) = int_0^1 H_+[g](Y) dY:

        delta_U(1) = (sin(pi*mu) / pi) * int v^{-mu} * ln(1 + 1/v) * G(1+v) dv.
    """
    n = len(v)
    w_xi = np.full(n, dxi)
    w_xi[0] *= 0.5
    w_xi[-1] *= 0.5
    row = np.zeros(n)
    pos = v > 0.0
    row[pos] = (
        (np.sin(np.pi * mu) / np.pi)
        * (v[pos] ** (-mu))
        * np.log1p(1.0 / v[pos])
        * dv_dxi[pos]
        * w_xi[pos]
    )
    return row


def L0_row(v: np.ndarray, dv_dxi: np.ndarray, dxi: float, mu: float) -> np.ndarray:
    """Row vector l0 such that l0 @ G = H_+[g](0) = (sin(pi*mu)/pi) * int v^{-mu}/(1+v) * G dv."""
    n = len(v)
    w_xi = np.full(n, dxi)
    w_xi[0] *= 0.5
    w_xi[-1] *= 0.5
    row = np.zeros(n)
    pos = v > 0.0
    row[pos] = (
        (np.sin(np.pi * mu) / np.pi)
        * (v[pos] ** (-mu))
        / (1.0 + v[pos])
        * dv_dxi[pos]
        * w_xi[pos]
    )
    return row


def build_operator_z(
    a: float,
    N: int = 1000,
    zmin: float = -14.0,
    zmax: float = 8.0,
    beta: float = 0.0,
    gauge: str = "a",
    edge_bc=None,
    edge_tail: bool = True,
    hilbert_grid_support=None,
    advection_scheme: str = "central",
) -> dict:
    """Assemble the desingularized linearized operator M_a on uniform z = ln(X-1) in [zmin, zmax].

    Parameters
    ----------
    a : float (< 0)
        gCLM advection parameter.
    N : int
        Number of grid points in z.
    zmin, zmax : float
        Log-coordinate domain bounds (v = X - 1 = e^z).
    beta : float
        Far-field exponential weight similarity W = diag(e^{beta * z}) (beta < 0 shifts
        transport essential spectrum left by (1-a)*beta).
    gauge : {'a', 'b'}
        'a': c_omega = -mu * c_l (kappa == 0, Hou-Luo Theta-balance analogue).
        'b': H(Omega)(0) == 1 (arXiv:2603.25104 Case 2.2 normalization).
    edge_bc : bool or None
        Impose G(z_min) = 0 by dropping node 0. If None (default): True for a >= -1
        (z_min is inflow), False for a < -1 (z_min is outflow!).
    edge_tail : bool
        Include exact sub-grid integral int_0^{v_0} H_+[g] dy in Q_op.
    hilbert_grid_support : bool or None
        Diagonal convention for hilbert_support_nonuniform (None follows edge_tail).
    advection_scheme : {'central', 'upwind2'}
        Finite-difference stencil for transport term.
    """
    a = float(a)
    mu = mu_of_a(a)
    c_l, c_w, _ = scaling_exponents(a)

    if edge_bc is None:
        edge_bc = bool(a >= -1.0)
    if hilbert_grid_support is None:
        hilbert_grid_support = bool(edge_tail)

    z = np.linspace(zmin, zmax, N)
    dz = z[1] - z[0]
    v = np.exp(z)
    X = 1.0 + v
    v_right = v[-1] + 0.5 * (v[-1] - v[-2])

    # Background advection speed in z: A_bar(1+v)/v = (1-a) * (1 - cos(pi*mu) * e^{-mu z})
    c_adv = (1.0 - a) * (1.0 - np.cos(np.pi * mu) * np.exp(-mu * z))

    if advection_scheme == "central":
        Dz = _d1_central(N, dz)
    elif advection_scheme == "upwind2":
        Dz = _d1_upwind2(N, dz, c_adv)
    else:
        raise ValueError(f"Unknown advection_scheme={advection_scheme!r}")

    # Hilbert transform acting on G via g = -sin(pi*mu) * v^{-mu} * G
    Hv = hilbert_support_nonuniform(v, v, dz, v_left=0.0, v_right=v_right,
                                    grid_support=hilbert_grid_support)
    g_scale = -np.sin(np.pi * mu) * (v ** (-mu))
    Hg_op = Hv * g_scale[None, :]  # (Hg_op @ G)_i = H_+[g](1 + v_i)

    # Exponential-averaging operator EA: (EA @ f)_i = int_{z_0}^{z_i} e^{z'-z_i} f(z') dz'
    EA = np.zeros((N, N))
    for i in range(1, N):
        w_ea = np.exp(z[: i + 1] - z[i])
        tw = np.full(i + 1, dz)
        tw[0] *= 0.5
        tw[i] *= 0.5
        EA[i, : i + 1] = w_ea * tw

    Q_op = EA @ Hg_op  # int_{v_0}^{v_i} H_+[g] dy / v_i
    t_row = edge_tail_row_z(z, mu) if edge_tail else np.zeros(N)
    if edge_tail:
        Q_op = Q_op + np.outer(np.exp(-z), t_row)

    # Transport + nonlocal bulk operator (Gauge a)
    transport = -c_adv[:, None] * Dz
    nonlocal_op = Hg_op + (a * mu) * Q_op
    M_bulk = transport + nonlocal_op

    # Linear functionals delta_U(1) and L_0[g] = H_+[g](0)
    u1_row = delta_U1_row(v, v, dz, mu)
    l0_row = L0_row(v, v, dz, mu)

    if gauge == "a":
        M = M_bulk.copy()
    elif gauge == "b":
        # kappa[G] = -l0_row @ (M_bulk @ G), so M = (I - 1 * l0_row) @ M_bulk
        l0_norm = float(np.sum(l0_row))  # ~ 1.0 (exact = 1 on infinite domain)
        M = M_bulk - np.outer(np.ones(N), (l0_row / l0_norm) @ M_bulk)
    else:
        raise ValueError(f"Unknown gauge={gauge!r}")

    # Store full-grid un-weighted operators for nonlinear RHS verification
    full_ops = dict(
        a=a, mu=mu, c_l=c_l, c_w=c_w, z=z, v=v, X=X, dz=dz,
        c_adv=c_adv, Dz=Dz, Hg_op=Hg_op, Q_op=Q_op, t_row=t_row,
        u1_row=u1_row, l0_row=l0_row, M_bulk=M_bulk, M_full=M,
        gauge=gauge, edge_bc=edge_bc,
    )

    if beta != 0.0:
        w_vec = np.exp(beta * z)
        wi_vec = np.exp(-beta * z)
        M = (w_vec[:, None] * M) * wi_vec[None, :]

    if edge_bc:
        M = M[1:, 1:]
        z_out = z[1:]
        v_out = v[1:]
        X_out = X[1:]
    else:
        z_out = z
        v_out = v
        X_out = X

    return dict(
        M=M,
        z=z_out,
        v=v_out,
        X=X_out,
        dz=dz,
        a=a,
        mu=mu,
        beta=beta,
        gauge=gauge,
        edge_bc=edge_bc,
        advection_scheme=advection_scheme,
        full=full_ops,
    )


def nonlinear_rhs_z(G: np.ndarray, ops: dict) -> np.ndarray:
    """Evaluate the RAW nonlinear dynamic-rescaling PDE RHS divided by Omega_bar_a:

        Omega_tau = -(c_l * X + a * U) * Omega_X + (c_omega + U_X) * Omega,
        G_tau     = Omega_tau / Omega_bar_a,

    for Omega = Omega_bar_a * (1 + G) on the full z-grid in `ops['full']`,
    using the exact same discrete operators (Dz, Hg_op, Q_op, u1_row, l0_row).
    Used for finite-difference certification of M_a (max|dG_lin - dG_fd| < 1e-10).
    """
    f = ops["full"]
    a = f["a"]
    mu = f["mu"]
    c_l_bar, c_w_bar = f["c_l"], f["c_w"]
    X = f["X"]
    v = f["v"]
    G = np.asarray(G, dtype=float)

    Om_b = omega_bar(X, a)
    OmX_b = omega_bar_X(X, a)
    U_b = U_bar(X, a)
    Ux_b = Ux_bar(X, a)

    # Perturbation fields from discrete operators
    Hg = f["Hg_op"] @ G
    # J(1+v) = int_1^{1+v} H_+[g] dy = v * (Q_op @ G)
    J = v * (f["Q_op"] @ G)
    dU1 = float(f["u1_row"] @ G)
    dU = dU1 + J  # since delta_U(X) = delta_U(1) + int_1^X H_+[g] dY

    # Full fields
    Om = Om_b * (1.0 + G)
    OmX = OmX_b * (1.0 + G) + Om_b * ((f["Dz"] @ G) / v)
    U = U_b + dU
    Ux = Ux_b + Hg

    # Edge-pinning normalization: c_l + a * U(1) = 0 -> c_l = c_l_bar - a * dU1
    c_l = c_l_bar - a * dU1

    if f["gauge"] == "a":
        c_w = -mu * c_l
    elif f["gauge"] == "b":
        # Choose c_w so that d/dtau H(Omega)(0) = 0, i.e., l0_row @ (Omega_tau / Omega_bar) = 0
        # Wait: Omega_tau = -(c_l X + a U)*OmX + (c_w + Ux)*Om
        # Compute c_w = 0 part first, then solve scalar linear eq for c_w!
        rhs_no_cw = -(c_l * X + a * U) * OmX + Ux * Om
        num = float(f["l0_row"] @ (rhs_no_cw / Om_b))
        den = float(f["l0_row"] @ (1.0 + G))
        c_w = -num / den
    else:
        raise ValueError(f"Unknown gauge={f['gauge']!r}")

    Om_tau = -(c_l * X + a * U) * OmX + (c_w + Ux) * Om
    return Om_tau / Om_b


def build_operator_zeta(
    a: float,
    N: int = 800,
    zetamax: float = 10.0,
    beta: float = 0.0,
    gauge: str = "a",
    edge_bc: bool = True,
    advection_scheme: str = "central",
) -> dict:
    """Assemble M_a for -1 < a < 0 in the EXACT characteristic coordinate zeta in [0, zetamax]:

        w = v^mu = kappa_a * (e^zeta - 1),    kappa_a = -cos(pi*mu) > 0,
        v(zeta) = (kappa_a * (e^zeta - 1))^{1/mu}.

    In this coordinate:
      * The physical singularity edge X = 1 is EXACTLY at zeta = 0 (no z_min cutoff!).
      * Background transport is CONSTANT UNIT SPEED: -A_bar_a * d_v = -d_zeta.
      * Inflow BC is G(0) = 0 (dropped DOF), so g = -sin(pi*mu)*v^{-mu}*G vanishes or is
        regularized cleanly at zeta = 0, and dv/dzeta = (1/mu)*kappa_a^{1/mu}*(e^zeta-1)^{1/mu - 1}*e^zeta
        vanishes at zeta = 0 since 1/mu = 1 - a > 1.
    """
    a = float(a)
    if not (-1.0 < a < 0.0):
        raise ValueError(f"zeta-coordinate requires -1 < a < 0 (cos(pi*mu) < 0); got a={a}")

    mu = mu_of_a(a)
    kappa = -np.cos(np.pi * mu)  # > 0
    inv_mu = 1.0 - a  # = 1/mu in (1, 2)

    zeta = np.linspace(0.0, zetamax, N)
    dzeta = zeta[1] - zeta[0]
    em1 = np.expm1(zeta)  # e^zeta - 1 >= 0
    w = kappa * em1
    v = w ** inv_mu
    X = 1.0 + v
    dv_dzeta = inv_mu * (kappa ** inv_mu) * (em1 ** (inv_mu - 1.0)) * np.exp(zeta)
    v_right = v[-1] + 0.5 * (v[-1] - v[-2])

    if advection_scheme == "central":
        Dzeta = _d1_central(N, dzeta)
    elif advection_scheme == "upwind2":
        Dzeta = _d1_upwind2(N, dzeta, np.ones(N))
    else:
        raise ValueError(f"Unknown advection_scheme={advection_scheme!r}")

    # Build Hilbert operator on interior nodes j = 1..N-1 (node 0 has v=0, G(0)=0, dv_dzeta=0)
    # Wait: include all N nodes with zero column 0 (since dv_dzeta[0] = 0 and G[0] = 0)
    Hv = np.zeros((N, N))
    # For i >= 1, j >= 1:
    v_int = v[1:]
    dv_int = dv_dzeta[1:]
    n_int = N - 1
    w_xi = np.full(n_int, dzeta)
    w_xi[-1] *= 0.5  # left endpoint is zeta_0=0 where integrand is 0 (full trapezoid weight dzeta at zeta_1!)
    w_v = dv_int * w_xi

    d = v_int[:, None] - v_int[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        K = w_v[None, :] / d
    K[np.diag_indices(n_int)] = 0.0
    # Exact p.v. of constant 1 over [0, v_right]: ln(v_i / (v_right - v_i))
    pv = np.log(v_int / (v_right - v_int))
    H_int = K / np.pi
    H_int[np.diag_indices(n_int)] = (pv - K.sum(axis=1)) / np.pi

    g_scale_int = -np.sin(np.pi * mu) * (v_int ** (-mu))
    Hg_int = H_int * g_scale_int[None, :]  # (n_int, n_int)

    # Q_op on interior nodes: Q_+[g](1 + v_i) = (1 / v_i) * int_0^{zeta_i} H_+[g] * dv_dzeta dzeta'
    # Trapezoid in zeta starting from zeta_0 = 0 where dv_dzeta(0) = 0 -> exact!
    Q_int = np.zeros((n_int, n_int))
    # Cumulative trapezoid matrix C on nodes 1..N-1 assuming integrand = 0 at node 0:
    # int_0^{zeta_1} f dzeta = 0.5*dzeta * f_1;  int_0^{zeta_i} = int_0^{zeta_{i-1}} + 0.5*dzeta*(f_{i-1}+f_i)
    C = np.zeros((n_int, n_int))
    for i in range(n_int):
        if i == 0:
            C[0, 0] = 0.5 * dzeta
        else:
            C[i] = C[i - 1]
            C[i, i - 1] += 0.5 * dzeta
            C[i, i] += 0.5 * dzeta
    # Ratio dv_dzeta[j] / v[i] (multiply on right by dv_int, on left by 1/v_int)
    Q_int = (C * (dv_int[None, :] / v_int[:, None])) @ Hg_int

    # Transport -d_zeta on interior nodes (with G(0) = 0 -> drop row/col 0)
    if not edge_bc:
        raise ValueError("zeta-coordinate always enforces physical inflow BC G(0)=0.")
    transport = -Dzeta[1:, 1:]
    M_bulk = transport + Hg_int + (a * mu) * Q_int

    if gauge == "a":
        M = M_bulk
    elif gauge == "b":
        # L0 functional on interior nodes
        l0_int = (
            (np.sin(np.pi * mu) / np.pi)
            * (v_int ** (-mu))
            / (1.0 + v_int)
            * w_v
        )
        l0_norm = float(np.sum(l0_int))
        M = M_bulk - np.outer(np.ones(n_int), (l0_int / l0_norm) @ M_bulk)
    else:
        raise ValueError(f"Unknown gauge={gauge!r}")

    if beta != 0.0:
        w_vec = np.exp(beta * zeta[1:])
        wi_vec = np.exp(-beta * zeta[1:])
        M = (w_vec[:, None] * M) * wi_vec[None, :]

    return dict(
        M=M,
        zeta=zeta[1:],
        v=v_int,
        X=X[1:],
        dzeta=dzeta,
        a=a,
        mu=mu,
        beta=beta,
        gauge=gauge,
        advection_scheme=advection_scheme,
    )


def spectrum(M: np.ndarray):
    """Compute eigenvalues and right eigenvectors sorted by decreasing real part."""
    vals, vecs = np.linalg.eig(M)
    idx = np.argsort(-vals.real)
    return vals[idx], vecs[:, idx]


def build_operator_branch2_z(
    a: float,
    N: int = 500,
    zmin: float = -14.0,
    zmax: float = 12.0,
    beta: float = 0.0,
    gauge: str = "a",
    edge_bc: bool = True,
    advection_scheme: str = "upwind2",
) -> dict:
    """Assemble the linearized operator M_{II,a} around Branch II (universal beta = 1/2):

        Omega_{II,a}(X) = (1 / a) * (X - 1)^{-1/2} * 1_{X > 1},   c_l = 2,   c_omega = -1.

    In z = ln(X - 1), with g = (1 / a) * v^{-1/2} * G:
        M_{II,a}^{bulk} G = -2 d_z G + (1 / a) H_+[v^{-1/2} G] + (1 / 2) Q_+[v^{-1/2} G].

    Reduces to M_{-1}^{bulk} at a = -1, and to 1D Hou-Luo (-2 d_z + (1/2) Q_+) as a -> -inf.
    """
    a = float(a)
    if a >= 0.0:
        raise ValueError(f"Expected a < 0; got a={a}")

    z = np.linspace(zmin, zmax, N)
    dz = z[1] - z[0]
    v = np.exp(z)
    X = 1.0 + v
    v_right = v[-1] + 0.5 * (v[-1] - v[-2])

    if advection_scheme == "central":
        Dz = _d1_central(N, dz)
    elif advection_scheme == "upwind2":
        Dz = _d1_upwind2(N, dz, np.full(N, 2.0))
    else:
        raise ValueError(f"Unknown advection_scheme={advection_scheme!r}")

    Hv = hilbert_support_nonuniform(v, v, dz, v_left=0.0, v_right=v_right, grid_support=True)
    # g = (1 / a) * v^{-1/2} * G
    Hg_op = Hv * ((1.0 / a) * (v ** -0.5))[None, :]

    EA = np.zeros((N, N))
    for i in range(1, N):
        w_ea = np.exp(z[: i + 1] - z[i])
        tw = np.full(i + 1, dz)
        tw[0] *= 0.5
        tw[i] *= 0.5
        EA[i, : i + 1] = w_ea * tw

    # edge_tail_row_z(z, 0.5) is built for g = -v^{-1/2}*G; rescale by (-1/a) for g = (1/a)*v^{-1/2}*G
    t_row = (-1.0 / a) * edge_tail_row_z(z, 0.5)
    Q_op = EA @ Hg_op + np.outer(np.exp(-z), t_row)

    M_bulk = -2.0 * Dz + Hg_op + (0.5 * a) * Q_op

    w_z = np.full(N, dz)
    w_z[0] *= 0.5
    w_z[-1] *= 0.5
    l0_row = (1.0 / np.pi) * (v ** 0.5) / (1.0 + v) * w_z  # H_+[g](0) / H_+[Om_II](0)

    if gauge == "a":
        M = M_bulk.copy()
    elif gauge == "b":
        M = M_bulk - np.outer(np.ones(N), (l0_row / np.sum(l0_row)) @ M_bulk)
    else:
        raise ValueError(f"Unknown gauge={gauge!r}")

    if beta != 0.0:
        w_vec = np.exp(beta * z)
        wi_vec = np.exp(-beta * z)
        M = (w_vec[:, None] * M) * wi_vec[None, :]

    if edge_bc:
        M = M[1:, 1:]
        z_out, v_out, X_out = z[1:], v[1:], X[1:]
    else:
        z_out, v_out, X_out = z, v, X

    return dict(
        M=M, z=z_out, v=v_out, X=X_out, dz=dz, a=a, beta=beta,
        gauge=gauge, edge_bc=edge_bc, advection_scheme=advection_scheme,
    )
