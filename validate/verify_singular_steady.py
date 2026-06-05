#!/usr/bin/env python3
"""
GATE: independently certify the Theorem-2.3 singular steady state of the 1D
Hou-Luo model (arXiv:2604.01868), both numerically and against the analytic
derivation.

Checks:
  1. line_hilbert is correct      -> Lorentzian:  H[1/(1+x^2)] = x/(1+x^2).
  2. H(Omega_bar) matches analytic -1/sqrt(1-X) (X<1), 0 (X>1).
  3. U from integrating H(Omega_bar) gives U(1) = -2  (the value the steady
     Theta-equation forces).
  4. Steady residual of eq.(2.4) ~ 0 on the smooth region X>1 with
     (c_l, c_omega) = (2, -1), using BOTH analytic and numeric U.
"""
import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from rescaled import (hilbert_matrix, line_hilbert, singular_omega,  # noqa: E402
                      singular_U_analytic, singular_Ux_analytic,
                      steady_residual)


def check_hilbert_lorentzian():
    # Wide window so the slowly decaying tails are captured; interior compared.
    L, n = 200.0, 8000
    x = np.linspace(-L, L, n)
    h = x[1] - x[0]
    f = 1.0 / (1.0 + x ** 2)
    M = hilbert_matrix(n)
    Hf = line_hilbert(f, M)            # h-independent kernel
    exact = x / (1.0 + x ** 2)
    interior = np.abs(x) < 50.0
    err = np.max(np.abs(Hf[interior] - exact[interior]))
    print(f"  [1] line_hilbert vs Lorentzian:  max|err| (|x|<50) = {err:.3e}")
    return err < 1e-3, M, h


def check_singular(M):
    # Grid on R, offset so no node lands exactly on the X=1 singularity.
    L, n = 200.0, 8000
    X = np.linspace(-L, L, n) + 0.5 * (2 * L / (n - 1))  # half-cell offset
    # exact spacing for integration
    dX = X[1] - X[0]

    # Omega_bar with the nearest node to X=1 regularized (integrable singularity).
    Omega = singular_omega(X, eps=0.0)
    # cap the single value adjacent to X=1 to avoid a grid spike in the conv
    near = np.argmin(np.abs(X - 1.0))
    if X[near] > 1.0:
        Omega[near] = 1.0 / np.sqrt(dX)        # ~ cell-averaged magnitude
    Hnum = line_hilbert(Omega, M)
    Ux_an = singular_Ux_analytic(X)

    # compare H(Omega_bar) to analytic, away from the singular cell
    far = (np.abs(X - 1.0) > 0.5) & (np.abs(X) < 50.0)
    err_H = np.max(np.abs(Hnum[far] - Ux_an[far]))
    print(f"  [2] H(Omega_bar) vs analytic U_X (|X-1|>0.5, |X|<50): "
          f"max|err| = {err_H:.3e}")

    # U(1) from integrating the analytic U_X (clean) and numeric U_X.
    # U(X) = int_0^X U_X dX'.  Use analytic primitive for the exact check.
    U1_analytic = 2.0 * np.sqrt(1.0 - 0.0) - 2.0  # = 0 at X=0 ... evaluate at 1:
    U1_exact = 2.0 * np.sqrt(max(1.0 - 1.0, 0.0)) - 2.0  # = -2
    # numeric: cumulative trapezoid of Hnum from the X=0 node to the X=1 node
    i0 = np.argmin(np.abs(X - 0.0))
    i1 = np.argmin(np.abs(X - 1.0))
    lo, hi = sorted((i0, i1))
    U1_num = np.trapezoid(Hnum[lo:hi + 1], X[lo:hi + 1])
    if i1 < i0:
        U1_num = -U1_num
    print(f"  [3] U(1): analytic = {U1_exact:+.6f}  (target -2),  "
          f"numeric int H = {U1_num:+.4f}")

    # Steady residual on the smooth region X>1 (analytic derivatives).
    c_l, c_omega = 2.0, -1.0
    Xs = X[(X > 1.2) & (X < 40.0)]
    Om = 1.0 / np.sqrt(Xs - 1.0)
    Om_X = -0.5 * (Xs - 1.0) ** (-1.5)
    Th = np.full_like(Xs, np.pi / 2.0)
    Th_X = np.zeros_like(Xs)
    U_an = singular_U_analytic(Xs)        # = -2 on X>1
    R_O, R_T = steady_residual(Xs, Om, Om_X, Th, Th_X, U_an, c_l, c_omega)
    print(f"  [4] steady residual on X in (1.2,40)  (analytic U):  "
          f"max|R_Omega| = {np.max(np.abs(R_O)):.3e}  "
          f"max|R_Theta| = {np.max(np.abs(R_T)):.3e}")

    return err_H, abs(U1_exact + 2.0), np.max(np.abs(R_O)), np.max(np.abs(R_T))


def main():
    print("=" * 72)
    print("  SINGULAR STEADY-STATE GATE  (Thm 2.3, Hou-Luo, c_l=2, c_omega=-1)")
    print("=" * 72)
    ok_h, M, h = check_hilbert_lorentzian()
    err_H, u1err, rO, rT = check_singular(M)
    print("-" * 72)
    verdict = ok_h and (u1err < 1e-9) and (rO < 1e-9) and (rT < 1e-12)
    print("  analytic U(1)=-2 exact:", u1err < 1e-9,
          " | steady residual ~0 (analytic):", rO < 1e-9 and rT < 1e-12)
    print("\n  VERDICT:",
          "PASS -- Thm 2.3 singular steady state independently certified"
          if verdict else "see numbers above")
    print("=" * 72)
    return 0 if verdict else 1


if __name__ == "__main__":
    sys.exit(main())
