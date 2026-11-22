#!/usr/bin/env python3
"""Probe: the constrained linearization is conjugate to damped transport.

Claim under test (docs/findings_session_2026_11_21_w3_conjugation.md):

  With N(s) = cos(pi s) Gamma(2s+2) / pi and c(s) = -tan(pi s)/(2s+1),

      N(s+1/2) c(s+1/2) = -N(s),

  so  M = -2 d_z + K  equals  N (-2 d_z - e^{-z/2}) N^{-1}  on L^2_sigma, 0 < sigma < 1/2,
  and the evolution is explicit.  For k0 = u^a e^{-b u}, u = e^{-z/2}:

      Ghat(s, tau) = N(s) e^{-2 s tau} 2 Gamma(2s+a) (b + x)^{-2s-a},   x = 1 - e^{-tau}.

Checks (exit 0 iff all pass):
  WC-1  the symbol identity, at complex s.
  WC-2  the closed form satisfies the ladder equation
        d/dtau Ghat(s) = -2 s Ghat(s) + c(s+1/2) Ghat(s+1/2)   (finite difference in tau).
  WC-3  the closed form agrees with expm(tau * M_GG) built from the repo's real-space
        discretization (src/stability_edge_z.py; p.v. Hilbert quadrature, no Mellin input),
        and the error falls when N doubles.  Transport alone is far off, as a contrast.

This probe does not test boundedness of the semigroup on L^2_sigma; that is an analytic
argument (the findings doc, section 2) plus the reviewer's operator-norm probe (section 4).
NumPy/SciPy only; about 30 s.
"""
import os
import sys

import numpy as np
from scipy.linalg import expm
from scipy.special import loggamma

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from stability_edge_z import build_blocks_z  # noqa: E402


def N_sym(s):
    return np.cos(np.pi * s) * np.exp(loggamma(2 * s + 2)) / np.pi


def c_sym(s):
    return -np.tan(np.pi * s) / (2 * s + 1)


A_EXP, B_PAR, SIGMA = 1.0, 1.0, 0.3


def G_hat(s, tau, a=A_EXP, b=B_PAR):
    x = 1.0 - np.exp(-tau)
    return N_sym(s) * np.exp(-2 * s * tau) * 2 * np.exp(loggamma(2 * s + a)) * (b + x) ** (-2 * s - a)


_T = np.linspace(-40.0, 40.0, 16001)
_DT = _T[1] - _T[0]
_SL = SIGMA + 1j * _T


def G_exact(z, tau):
    return (np.exp(np.outer(z, _SL)) @ G_hat(_SL, tau)).real * _DT / (2 * np.pi)


def main():
    ok = True

    s = np.array([0.3 + 0.7j, 0.1 - 3j, 0.45 + 12j, -0.2 + 0.4j, 0.25])
    e1 = np.max(np.abs(N_sym(s + 0.5) * c_sym(s + 0.5) + N_sym(s)) / np.abs(N_sym(s)))
    p1 = e1 < 1e-12
    ok &= p1
    print(f"WC-1 symbol identity             max rel err {e1:.2e}   {'PASS' if p1 else 'FAIL'}")

    h, tau = 1e-5, 0.7
    s = np.array([0.3 + 0.5j, 0.2 - 2j, 0.4 + 6j])
    lhs = (G_hat(s, tau + h) - G_hat(s, tau - h)) / (2 * h)
    rhs = -2 * s * G_hat(s, tau) + c_sym(s + 0.5) * G_hat(s + 0.5, tau)
    e2 = np.max(np.abs(lhs - rhs) / np.abs(rhs))
    p2 = e2 < 1e-6
    ok &= p2
    print(f"WC-2 ladder equation             max rel err {e2:.2e}   {'PASS' if p2 else 'FAIL'}")

    errs = {}
    for n in (1500, 3000):
        blk = build_blocks_z(N=n, zmin=-14.0, zmax=12.0, edge_bc=True)
        z, m = blk["z"], blk["M_GG"]
        g0 = G_exact(z, 0.0)
        win = (z > -8) & (z < 5)
        worst, worst_tr = 0.0, np.inf
        for tau in (0.25, 0.5, 1.0):
            gd = expm(tau * m) @ g0
            worst = max(worst, np.abs(gd - G_exact(z, tau))[win].max())
            worst_tr = min(worst_tr, np.abs(gd - G_exact(z - 2 * tau, 0.0))[win].max())
        errs[n] = worst
        print(f"     N={n}: max|G0|={np.abs(g0).max():.4f}  err vs closed form {worst:.2e}"
              f"  (vs transport alone {worst_tr:.2e})")
    p3 = errs[1500] < 1e-3 and errs[3000] < 0.7 * errs[1500]
    ok &= p3
    print(f"WC-3 real-space discretization   {errs[1500]:.2e} -> {errs[3000]:.2e}   {'PASS' if p3 else 'FAIL'}")

    print("ALL PASS" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
