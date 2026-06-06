#!/usr/bin/env python3
"""
CO-EVOLVING-base dynamic-rescaling march (THE NEXT CHIP #2), STARTED.

Goal: let c_l, c_omega adjust dynamically by the normalization and confirm a
perturbed profile RELAXES to (Omega_bar, c_l=2, c_omega=-1) with no exponent drift.

This script (a) validates the dynamic-normalization machinery in the VERIFIED split
form (closed-form base velocity Ubar + bounded perturbation Hilbert on G), so the
profile is an EXACT fixed point of the normalization (c_l=2, c_omega=-1 to machine
precision), and (b) PINPOINTS + QUANTIFIES the wall: in the fixed-Y frame the c_l
(dilation) gauge couples to the edge-singular mode  X*Omega_bar_X ~ (Y-1)^{-3/2}
(more singular than the profile), whose G-variable magnitude ~ e^{-z} DIVERGES at the
edge and grows with resolution -> the fixed-Y frame is INVALID; the moving edge-fixed
frame (xi-ODE, findings_edge_condition_derivation.md) is mandatory.

Split (verified) variables:  Omega = Omega_bar + dOmega,  dOmega = (Y-1)^{-1/2} G,
G = W-1 bounded;  Theta = Theta_bar + h.  Closed form  Ubar=-2 on Y>1, Ubar(1)=-2.
  dU(Y>1) = int_1^Y H(dOmega) dY' = J_op @ G        (verified, bounded)
  dU(1)   = int_0^1 H(dOmega) dY' = b1 @ G          (Hilbert at Y<1, regular)
  c_l = -U(1) = -(Ubar(1)+dU(1)) = 2 - dU(1)        (= 2 at profile, G=0)
  c_omega = [abar*Omega_Y - Theta_Y]/Omega |_Yref   (= -1 at profile)
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from stability_edge_z import build_blocks_z, hilbert_support_z  # noqa: E402


def build(N=1400, zmin=-12.0, zmax=12.0, Nm=600):
    B = build_blocks_z(N=N, zmin=zmin, zmax=zmax, beta=0.0, edge_bc=False)
    z, dz, Y = B["z"], B["dz"], B["Y"]
    ym1 = Y - 1.0
    Hg_op, J_op, Dz = B["Hg_op"], B["J_op"], B["Dz"]

    # b1 @ G = dU(1) = int_0^1 H(dOmega) dY' ,  dOmega = e^{-z/2} G  (Y'<1<Y, regular)
    Ym = np.linspace(0.0, 1.0, Nm)
    wYsupp = np.empty(N)
    wYsupp[1:-1] = 0.5 * (Y[2:] - Y[:-2])
    wYsupp[0] = 0.5 * (Y[1] - Y[0])
    wYsupp[-1] = 0.5 * (Y[-1] - Y[-2])
    Hm = (1.0 / np.pi) * (wYsupp[None, :] / (Ym[:, None] - Y[None, :]))   # H(.) at Y<1
    Hm_G = Hm @ np.diag(ym1 ** -0.5)                                      # acts on G
    # int_0^1 over Ym (trapezoid), U(0)=0
    tw = np.full(Nm, Ym[1] - Ym[0]); tw[0] *= 0.5; tw[-1] *= 0.5
    b1 = tw @ Hm_G
    return dict(z=z, dz=dz, Y=Y, ym1=ym1, Hg_op=Hg_op, J_op=J_op, Dz=Dz, b1=b1, N=N,
                zmin=zmin)


def normalization(P, G, h, iref):
    """Return c_l, c_omega from the dynamic normalization (split, exact at profile)."""
    z, Y, ym1, Dz = P["z"], P["Y"], P["ym1"], P["Dz"]
    dU1 = float(P["b1"] @ G)                       # dU(1)
    c_l = 2.0 - dU1                                # = -U(1)
    dU = P["J_op"] @ G                             # dU(Y) on Y>1
    U = -2.0 + dU                                  # Ubar=-2 on Y>1
    abar = U + c_l * Y
    # Omega = e^{-z/2}(1+G);  Omega_Y = e^{-3z/2}(G_z - (1+G)/2)
    Gz = Dz @ G
    Omega = np.exp(-z / 2) * (1.0 + G)
    Omega_Y = np.exp(-1.5 * z) * (Gz - 0.5 * (1.0 + G))
    Theta_Y = np.exp(-z) * (Dz @ h)                # Theta_bar const -> only h
    i = iref
    c_omega = float((abar[i] * Omega_Y[i] - Theta_Y[i]) / (Omega[i] + 1e-300))
    return c_l, c_omega, dU1, abar


def main():
    print("=" * 80)
    print("  CO-EVOLVING march STARTED -- dynamic normalization + the wall (NEXT CHIP #2)")
    print("=" * 80)

    P = build(N=1400, zmin=-12.0, zmax=12.0)
    z, Y, ym1 = P["z"], P["Y"], P["ym1"]
    iref = int(np.argmin(np.abs(Y - 3.0)))

    # ---- GATE 0a: validate the dU(1) functional against direct fine quadrature ----
    Gtest = np.exp(-((z - 1.5) / 1.2) ** 2)        # smooth bounded test G
    dU1_op = float(P["b1"] @ Gtest)
    # direct: dOmega = (Y-1)^{-1/2} Gtest on Y>1; H at a fine Y'<1 grid; integrate
    Yf = np.linspace(0.0, 1.0, 4000)
    wY = np.empty(P["N"]); wY[1:-1] = 0.5 * (Y[2:] - Y[:-2])
    wY[0] = 0.5 * (Y[1] - Y[0]); wY[-1] = 0.5 * (Y[-1] - Y[-2])
    dOm = ym1 ** -0.5 * Gtest
    Hf = (1 / np.pi) * (wY[None, :] / (Yf[:, None] - Y[None, :])) @ dOm
    dU1_direct = float(np.trapz(Hf, Yf)) if hasattr(np, "trapz") else float(np.trapezoid(Hf, Yf))
    print("\n  GATE 0a -- dU(1) functional vs direct quadrature (validate the new piece):")
    print(f"    b1@G = {dU1_op:+.5f}   direct int_0^1 H(dOmega) = {dU1_direct:+.5f}   "
          f"rel.err {abs(dU1_op-dU1_direct)/(abs(dU1_direct)+1e-30):.2e}")

    # ---- GATE 0b: normalization on the EXACT profile (G=0, h=0) -> (2,-1) exact ----
    G0 = np.zeros(P["N"]); h0 = np.zeros(P["N"])
    c_l, c_om, dU1, abar = normalization(P, G0, h0, iref)
    print("\n  GATE 0b -- normalization on the exact profile (split form):")
    print(f"    dU(1) = {dU1:+.2e}   c_l = {c_l:+.8f} (target 2)   "
          f"c_omega = {c_om:+.8f} (target -1)   abar(edge)={abar[0]:+.2e}")

    # ---- THE WALL: dilation gauge term is edge-singular in the fixed-Y frame ------
    # The c_l-correction forcing in eq.(2.4) is  -dc_l * X * Omega_X ; in the G-variable
    # (multiply by (Y-1)^{1/2}) the base part is  -dc_l * X * (Y-1)^{1/2} Omega_bar_X
    #   = -dc_l * (1+e^z) * e^{z/2} * (-1/2)(Y-1)^{-3/2} = (dc_l/2)(1+e^z) e^{-z}.
    # ~ (dc_l/2) e^{-z}  at the edge  -> DIVERGES as z->-inf (resolution-dependent).
    print("\n  THE WALL -- fixed-Y dilation gauge term  -dc_l * X * Omega_bar_X, in G:")
    print("    (perturb with a bulk bump so dc_l = -dU(1) != 0; report edge magnitude)")
    print("    zmin     N      dc_l       max|gauge term|   value at edge node")
    for (Nc, zc) in [(1000, -10.0), (1400, -12.0), (1800, -14.0), (2200, -16.0)]:
        Pc = build(N=Nc, zmin=zc, zmax=12.0)
        zc_arr, ym1c = Pc["z"], Pc["ym1"]
        Gp = 0.15 * np.exp(-((zc_arr - 1.0) / 1.0) ** 2)
        dc_l = -float(Pc["b1"] @ Gp)
        Xc = 1.0 + np.exp(zc_arr)
        gauge_term = -dc_l * Xc * np.exp(0.5 * zc_arr) * (-0.5) * ym1c ** -1.5   # in G-var
        print(f"    {zc:+.0f}    {Nc:<5d}  {dc_l:+.4f}    {np.max(np.abs(gauge_term)):.3e}"
              f"        {gauge_term[0]:+.3e}")

    print("\n" + "=" * 80)
    print("  READ:")
    print("  - GATE 0a/0b: the dynamic-normalization machinery is correct -- the profile")
    print("    is an EXACT fixed point (c_l=2, c_omega=-1 to machine precision) in the")
    print("    verified split form.  (The earlier raw-Hilbert version gave 2.24/-1.11: an")
    print("    edge under-resolution artifact from Hilberting the singular profile.)")
    print("  - THE WALL: in the FIXED-Y frame the c_l (dilation) gauge term scales as")
    print("    (dc_l/2) e^{-z} at the edge -- it DIVERGES as z_min->-inf and grows with")
    print("    resolution.  X*Omega_bar_X ~ (Y-1)^{-3/2} is more singular than the profile")
    print("    (= the dilation/edge symmetry generator, weight-excluded in the linear")
    print("    verdict).  So a co-evolving march CANNOT run in the fixed-Y frame; the")
    print("    MOVING edge-fixed frame (xi-ODE: xi_dot=(c_l-1)xi+V_g(1)+mu) is mandatory.")
    print("    That is the precise next step (equations already derived in")
    print("    findings_edge_condition_derivation.md; the LINEAR operator M_GG already")
    print("    uses it).  The co-evolving NONLINEAR march = nonlinear_shape_march + xi-ODE")
    print("    + explicit mu,nu normalization, in the moving frame.")
    print("=" * 80)


if __name__ == "__main__":
    main()
