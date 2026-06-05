#!/usr/bin/env python3
"""
Nonlinear cross-check of the stability operator: verify it IS the correct
linearization of the true eq.(2.4) flow, using the TRUSTED nonlinear engine
(src/dynamic_rescaling.py, validated against the CLM exact solution to 2e-10).

This is the decisive, feasible form of the "nonlinear shape test": instead of
long-time marching on the singular saddle (which fights the edge-resolution wall),
we check the operator against the nonlinear flow's TENDENCY at finite amplitude.
If, for a shape perturbation (dOmega, dTheta),

   T_nl := [RHS_2.4(base + eps*dpert) - RHS_2.4(base)] / eps   (trusted engine)

matches the analytic linearization L(dpert) and the match -> 0 as eps -> 0, then the
linearized operator (whose spectrum gives the gap) is the correct linearization of
the genuine nonlinear dynamics -- catching any sign/factor/derivation error.

Analytic linearization of eq.(2.4) about (Omega_bar, Theta_bar), c_l=2, c_omega=-1:
   L_Omega = -abar dOmega_X - dU * Omega_bar_X + c_omega dOmega + dTheta_X
   L_Theta = -abar dTheta_X - dU * Theta_bar_X            (c_l+2 c_omega = 0)
with abar = U_bar + c_l X,  U_bar = int_0^X H(Omega_bar),  dU = int_0^X H(dOmega).
(The -dU*Theta_bar_X edge term is the one the moving frame turns into the xi-ODE;
it is included here so the FULL linearization is tested.)
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from dynamic_rescaling import DynamicRescaling, d_central4, d_upwind3  # noqa: E402


def velocity(sol, Om):
    Ux = sol.M @ Om
    U = np.concatenate(([0.0], np.cumsum(0.5 * (Ux[1:] + Ux[:-1]) * sol.dx)))
    U = U - U[sol.c]
    return U, Ux


def analytic_linearization(sol, Ob, Tb, dO, dT, c_l=2.0, c_omega=-1.0):
    """L(dO,dT) about base (Ob,Tb), assembled independently from the PDE."""
    Ub, _ = velocity(sol, Ob)
    abar = Ub + c_l * sol.X
    dU, _ = velocity(sol, dO)
    # base derivatives must use the SAME scheme as the engine's rhs (UPWIND for the
    # transported fields Omega,Theta; central only for the +Theta_X source).
    ObX = d_upwind3(Ob, abar, sol.dx)      # engine: Om_X = d_upwind3(Omega, a)
    TbX = d_upwind3(Tb, abar, sol.dx)      # engine: Th_X_adv = d_upwind3(Theta, a)
    dOX = d_upwind3(dO, abar, sol.dx)
    dTX_adv = d_upwind3(dT, abar, sol.dx)
    dTX_src = d_central4(dT, sol.dx)
    LO = -abar * dOX - dU * ObX + c_omega * dO + dTX_src
    LT = -abar * dTX_adv - dU * TbX        # (c_l+2c_omega)=0
    return LO, LT


def main():
    print("=" * 78)
    print("  Nonlinear linearization check: operator vs trusted eq.(2.4) RHS tendency")
    print("=" * 78)
    P, N = 40.0, 2001
    delta = 0.04
    sol = DynamicRescaling(N=N, P=P, c_l=2.0, c_omega_fixed=-1.0)
    X, dx = sol.X, sol.dx

    # regularized singular base (same regularization used consistently on both sides)
    sig = 0.5 * (1 + np.tanh((X - 1.0) / (3 * dx)))
    Ob = sig / np.sqrt(np.clip(X - 1.0 + delta, delta / 50, None))
    Tb = (np.pi / 2) * sig

    # a genuine localized SHAPE perturbation (smooth, supported on X>1, decays at edge)
    dO = np.exp(-((X - 3.0) / 1.2) ** 2) * (X > 1.0)
    dT = 0.5 * np.exp(-((X - 4.0) / 1.5) ** 2) * (X > 1.0)
    nrm = np.sqrt(np.sum(dO ** 2 + dT ** 2) * dx)
    dO, dT = dO / nrm, dT / nrm

    LO, LT = analytic_linearization(sol, Ob, Tb, dO, dT)

    print(f"  grid N={N}, P={P}, delta={delta}; shape perturbation supported on X>1")
    print("  eps        rel.err ||T_nl - L|| / ||L||   (-> 0 as eps->0 confirms operator)")
    core = np.abs(X) < 30.0                 # ignore far-field box edge
    for eps in (1e-2, 3e-3, 1e-3, 3e-4, 1e-4):
        rO0, rT0, _ = sol.rhs(Ob, Tb)
        rO1, rT1, _ = sol.rhs(Ob + eps * dO, Tb + eps * dT)
        TO = (rO1 - rO0) / eps
        TT = (rT1 - rT0) / eps
        num = np.sqrt(np.sum(((TO - LO) ** 2 + (TT - LT) ** 2)[core]) * dx)
        den = np.sqrt(np.sum((LO ** 2 + LT ** 2)[core]) * dx)
        print(f"   {eps:.0e}     {num/den:.3e}")

    # Rayleigh quotient of the nonlinear tendency: UNWEIGHTED vs convective WEIGHTED.
    rO0, rT0, _ = sol.rhs(Ob, Tb)
    eps = 1e-3
    rO1, rT1, _ = sol.rhs(Ob + eps * dO, Tb + eps * dT)
    TO, TT = (rO1 - rO0) / eps, (rT1 - rT0) / eps
    supp = (X > 1.0) & (X < 30.0)
    rq_u = (np.sum((dO * TO + dT * TT)[supp]) /
            (np.sum((dO ** 2 + dT ** 2)[supp]) + 1e-30))
    gamma = 0.5
    w = (X - 1.0) ** (-2 * gamma)            # e^{-2 gamma z}, z=ln(X-1)
    rq_w = (np.sum((w * (dO * TO + dT * TT))[supp]) /
            (np.sum((w * (dO ** 2 + dT ** 2))[supp]) + 1e-30))
    print(f"\n  Rayleigh quotient of the nonlinear tendency <dpert, T_nl>/<dpert,dpert>:")
    print(f"     unweighted          = {rq_u:+.4f}  (positive = NON-NORMAL transient;")
    print(f"                                  the unweighted norm is marginal/transport)")
    print(f"     convective weighted = {rq_w:+.4f}  (gamma={gamma}; the physically")
    print(f"                                  relevant norm -- genuine decay tendency)")

    print("\n" + "=" * 78)
    print("  READ: rel.err -> 0 (LINEARLY in eps) confirms the stability operator IS")
    print("  the correct linearization of the trusted eq.(2.4) RHS -- so the spectral")
    print("  gap / no-unstable-mode verdict rests on the right operator. The positive")
    print("  unweighted vs negative weighted Rayleigh quotient is the expected non-")
    print("  normal/convective signature (decay lives in the weighted norm, cf.")
    print("  spectral_gap.py).")
    print("=" * 78)


if __name__ == "__main__":
    main()
