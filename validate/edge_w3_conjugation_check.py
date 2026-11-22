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
  WC-4  the closed formula for Re c(x+it), negative for 0 < x <= 1/2 and positive near t = 0
        for 1/2 < x < 1.
  WC-5  the energy-form identity
        Re<(M+2 sigma)G, G>_sigma = (1/2pi) int Re c(sigma+1/4+it) |Ghat(sigma+1/4+it)|^2 dt,
        left side on the line Re s = sigma, right side on Re s = sigma + 1/4; non-positive for
        sigma <= 1/4, positive for suitable data when sigma > 1/4.
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

    # WC-4: closed formula for Re c(x+it) and its sign on 0 < x <= 1/2.
    def re_c_formula(x, t):
        num = (2 * x + 1) * np.sin(2 * np.pi * x) + 2 * t * np.sinh(2 * np.pi * t)
        den = (np.cos(2 * np.pi * x) + np.cosh(2 * np.pi * t)) * np.abs(2 * (x + 1j * t) + 1) ** 2
        return -num / den

    tt = np.linspace(-6.0, 6.0, 2400)          # even count: t = 0 is not a node
    e4, sign_ok = 0.0, True
    for x in (0.05, 0.2, 0.35, 0.45, 0.5, 0.6, 0.7):
        direct = c_sym(x + 1j * tt).real
        e4 = max(e4, np.max(np.abs(direct - re_c_formula(x, tt)) / (1 + np.abs(direct))))
        if x <= 0.5:
            sign_ok &= bool(np.all(direct <= 1e-13))
        else:
            sign_ok &= bool(direct[np.argmin(np.abs(tt))] > 0)
    p4 = e4 < 1e-9 and sign_ok
    ok &= p4
    print(f"WC-4 Re c formula and sign       max err {e4:.2e}  sign {'ok' if sign_ok else 'BAD'}   {'PASS' if p4 else 'FAIL'}")

    # WC-5: energy-form identity
    #   Re <(M+2 sigma) G, G>_sigma = (1/2pi) int Re c(sigma+1/4+it) |Ghat(sigma+1/4+it)|^2 dt
    # with the left side computed on the line Re s = sigma from (MG)^ = N(s)(-2s khat(s) - khat(s+1/2)).
    def k_hat(s, ell, s0):
        return np.exp(ell * (s - s0) ** 2)

    t5 = np.linspace(-30.0, 30.0, 24000)
    d5 = t5[1] - t5[0]
    e5, signs = 0.0, []
    for sig, ell, s0 in ((0.1, 4.0, 0.3), (0.2, 20.0, 0.45), (0.25, 6.0, 0.5), (0.3, 20.0, 0.55), (0.4, 20.0, 0.65)):
        sl = sig + 1j * t5
        g = N_sym(sl) * k_hat(sl, ell, s0)
        mg = N_sym(sl) * (-2 * sl * k_hat(sl, ell, s0) - k_hat(sl + 0.5, ell, s0))
        lhs = (np.sum(mg * np.conj(g)).real + 2 * sig * np.sum(np.abs(g) ** 2)) * d5 / (2 * np.pi)
        w = sig + 0.25 + 1j * t5
        gw = N_sym(w) * k_hat(w, ell, s0)
        rhs = np.sum(c_sym(w).real * np.abs(gw) ** 2) * d5 / (2 * np.pi)
        nrm = np.sum(np.abs(g) ** 2) * d5 / (2 * np.pi)
        e5 = max(e5, abs(lhs - rhs) / max(abs(rhs), 1e-300))
        signs.append((sig, rhs / nrm))
    neg_ok = all(v <= 0 for sg, v in signs if sg <= 0.25)
    pos_seen = any(v > 0 for sg, v in signs if sg > 0.25)
    p5 = e5 < 1e-6 and neg_ok and pos_seen
    ok &= p5
    print("     form/||G||^2: " + ", ".join(f"sigma={sg}: {v:+.4f}" for sg, v in signs))
    print(f"WC-5 energy-form identity        max rel err {e5:.2e}   {'PASS' if p5 else 'FAIL'}")

    print("ALL PASS" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
