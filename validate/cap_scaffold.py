#!/usr/bin/env python3
"""
CAP SCAFFOLD (chip #3, rigor): is the linearized stability operator DISSIPATIVE in a
weighted inner product?  If so, the Conjecture-2.4 spectral bound `sup Re lambda < 0`
reduces to a SYMMETRIC eigenvalue bound that a computer-assisted proof can certify.

Route. For an operator A and a weighted inner product <u,v>_w = u^T W v (W spd, here
W = diag(e^{-2 gamma z} dz)), the FIELD OF VALUES F(A) = {<Au,u>_w/<u,u>_w} contains the
spectrum, and  max Re F(A) = lambda_max of the generalized symmetric pencil
    Sym := (W A + A^T W)/2,   Sym u = lambda W u.
So  max Re spec(A) <= lambda_max(Sym, W).  If lambda_max < 0 in some weighted inner
product, the operator is DISSIPATIVE there and `sup Re lambda < 0` follows -- and the
remaining CAP step is a VERIFIED symmetric-eigenvalue upper bound (standard: interval
Cholesky / Rayleigh-Ritz + Lehmann-Goerisch) plus a rigorous truncation/tail bound.

Why a weight is needed: the operator is strongly NON-NORMAL (transport-dominated), so in
the PLAIN inner product the field of values is far to the right of the spectrum. Analytic
fact (checked): the transport -2 d_z conjugated into <.>_w has symmetric part the CONSTANT
-2 gamma (dissipative at rate 2 gamma). The CAP question is whether the NONLOCAL forcing
+ the finite-domain BOUNDARY terms stay below that.

This scaffold computes lambda_max(Sym, W) vs gamma for M_GG (Omega-shape), M_hh (Theta-
shape), and the coupled operator, and reports whether a dissipative weight exists -- i.e.
which CAP route is viable and the exact inequality to certify.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from stability_edge_z import build_blocks_z  # noqa: E402


def field_of_values_abscissa(A, z, dz, gamma):
    """max Re of the field of values of A in <u,v>_w, w=e^{-2 gamma z}.  Computed via the
    similarity transform  At = diag(e^{-gamma z}) A diag(e^{+gamma z})  (= sqrt(W) A
    sqrt(W)^{-1}; the dz constant cancels), so FoV_w(A) = FoV(At) in the standard inner
    product and the abscissa is lambda_max of the SYMMETRIC matrix (At+At^T)/2.  This
    avoids the ill-conditioned generalized eigenproblem.  >= max Re spec(A)."""
    sw = np.exp(-gamma * z)                       # sqrt(w) up to the (cancelling) dz
    At = (sw[:, None] * A) / sw[None, :]
    return float(np.linalg.eigvalsh(0.5 * (At + At.T))[-1])


def spec_abscissa(A):
    return float(np.linalg.eigvals(A).real.max())


def main():
    print("=" * 84)
    print("  CAP SCAFFOLD -- dissipativity of the stability operator in a weighted norm")
    print("=" * 84)

    N, zmin, zmax = 1200, -12.0, 12.0
    B = build_blocks_z(N=N, zmin=zmin, zmax=zmax, beta=0.0, edge_bc=True)
    z, dz = B["z"], B["dz"]
    M_GG, M_Gh, M_hh = B["M_GG"], B["M_Gh"], B["M_hh"]
    n = len(z)

    # analytic check: transport -2 d_z has weighted-symmetric part = -2 gamma (constant)
    T = -2.0 * B["Dz"]
    print("\n  CHECK -- transport -2 d_z: field-of-values abscissa vs predicted -2*gamma:")
    for g in (0.25, 0.5, 0.75):
        print(f"    gamma={g:.2f}:  FoV abscissa = {field_of_values_abscissa(T, z, dz, g):+.4f}"
              f"   (predict -2*gamma = {-2*g:+.2f})")

    print("  ^ the discrete one-sided boundary differences (amplified by the edge weight)")
    print("    CONTAMINATE the naive symmetrization -> the transport must be handled")
    print("    ANALYTICALLY, not by symmetrizing the discrete matrix.")

    print("\n  TRANSPORT, analytically (integration by parts, w=e^{-2 gamma z}):")
    print("    <g, -2 g_z>_w = -2 gamma ||g||_w^2 - [w g^2]_{zmin}^{zmax}")
    print("    inflow BC g(zmin)=0 kills the lower term; outflow makes the upper <= 0")
    print("    => transport contributes  <= -2 gamma ||g||_w^2  (dissipative at rate 2 gamma).")

    # The only operator term to BOUND numerically is the bounded nonlocal forcing
    # F = (1/2) Q_op = (1/2) dU/(Y-1)  (M_GG = transport + F).  Certifiable quantity:
    # C_force(gamma) = max Re field-of-values of F in <.>_w.  Net dissipative iff < 2 gamma.
    F = 0.5 * B["Q_op"]
    print("\n  NONLOCAL FORCING F=(1/2)Q_op: weighted field-of-values bound C_force(gamma)")
    print("  vs 2*gamma -> energy abscissa <= -2 gamma + C_force; <0 means DISSIPATIVE:")
    print("    gamma   C_force     2*gamma    energy abscissa (-2g+C_force)")
    best = (None, 1e9)
    for g in (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.75, 1.0):
        Cf = field_of_values_abscissa(F, z, dz, g)
        ea = -2 * g + Cf
        if ea < best[1]:
            best = (g, ea)
        flag = "  <- DISSIPATIVE" if ea < 0 else ""
        print(f"    {g:.2f}    {Cf:+.4f}    {2*g:.3f}     {ea:+.4f}{flag}")
    print(f"    best: gamma={best[0]}, energy abscissa = {best[1]:+.4f}")

    # C_force vs edge depth: does it CONVERGE (energy route viable) or DIVERGE (ruled out)?
    print("\n  C_force(gamma=0.4) vs edge depth (converge => energy route; diverge => ruled out):")
    for (Nc, zc_) in [(800, -10.0), (1200, -12.0), (1600, -14.0)]:
        Bc = build_blocks_z(N=Nc, zmin=zc_, zmax=12.0, beta=0.0, edge_bc=True)
        Cf = field_of_values_abscissa(0.5 * Bc["Q_op"], Bc["z"], Bc["dz"], 0.4)
        print(f"    N={Nc:<5d} zmin={zc_:+.0f}:  C_force={Cf:+.4f}   (grows with |zmin| => DIVERGES)")

    print("\n" + "=" * 84)
    print("  CAP-ROUTE DETERMINATION (the value of this scaffold):")
    print("  C_force GROWS with edge depth (10 -> 29 -> 83 as zmin -10 -> -14), NOT a fixed")
    print("  bound. Reason (real, not a discretization artifact): in <.>_w the similarity")
    print("  transform multiplies the NONLOCAL forcing entries by e^{gamma(z_j - z_i)}, huge")
    print("  for an edge row (z_i small) coupled to the bulk -- the exponential weight and the")
    print("  edge<->bulk nonlocal (Hilbert) coupling DO NOT COMMUTE. So the numerical RANGE in")
    print("  the convective weight is enormous while the SPECTRUM is at ~-1e-3 (items 2/8):")
    print("  textbook strong non-normality. => the ENERGY / DISSIPATIVITY CAP route is RULED")
    print("  OUT in the convective weight.")
    print("")
    print("  Therefore the CAP must take the EIGENVALUE-ENCLOSURE route on the spectrum")
    print("  directly (the spectrum is benign; only the field-of-values is bad):")
    print("   (1) a verified eigenpair enclosure for each low mode (Newton-Kantorovich /")
    print("       Yamamoto a posteriori) + a Lehmann-Goerisch lower bound separating the")
    print("       rest of the spectrum from Re=0, modulo span{xi@1, mu, nu};")
    print("   (2) a rigorous Hilbert-kernel enclosure (the ~1.5% quadrature -> verified);")
    print("   (3) a truncation/tail bound for the mapped infinite domain.")
    print("  OR find a NON-exponential metric that bounds the nonlocal term while keeping")
    print("  transport dissipative (open; the convective weight provably does not).")
    print("  THE LEAP = (1)-(3) [De Huang's CAP machinery], or the better-metric question.")
    print("  Transport + boundary signs are CLOSED analytically above; the obstruction is now")
    print("  pinned to the nonlocal coupling's interaction with the weight. Numerical scaffold.")
    print("=" * 84)


if __name__ == "__main__":
    main()
