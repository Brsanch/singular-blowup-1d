#!/usr/bin/env python3
"""
OBS-07 gate: depth / weak-coupling scan of the exact ladder resolvent in the convective
weight L2_sigma (docs/desk/convective_weight_spectrum.md Section 7).

For F(z) = e^{-cosh(z - z0)} translated edge-ward, computes
    gain(lam, z0) = ||G||_sigma / ||F||_sigma,   G = (lam - M)^{-1} F
via the defect-closed 4-Gamma ladder (LEM-09; machine-verified against the continuum
operator in validate/edge_convective_structure.py). By the eps-similarity this scans the
weak-coupling family eps = e^{z0/2}: bounded-in-z0 dressing supports
sigma(M|_sigma) inside {Re lam <= -2 sigma} with resolvent bounds (the essential-spectrum
half of LEM-06); growth like e^{theta|z0|} would refute the convective-weight gap.

Numerical notes: the ladder terms alternate and reach ~e^{2 e^{|z0|/4}} before cancelling —
that is WHY termwise analytic bounds fail (desk doc Section 8) and why this runs at high dps.
Everything is Gamma/Mellin arithmetic on the line Re s = sigma: NO discrete Hilbert, immune
to the repo's known edge-bias trap.

Default = quick mode (z0 down to -10, two lambdas, dps 60, ~minutes).
  --deep : the full grid of the desk doc (z0 to -14, six lambdas, dps 90).
PASS criteria (structural, not stability-asserting): the line-equation self-check of the
constructed Ghat and finiteness of every gain.
"""
import sys
import time

import mpmath as mp

DEEP = "--deep" in sys.argv
DPS = 90 if DEEP else 60
mp.mp.dps = DPS
SIGMA = mp.mpf('0.4')


def c_sym(s):
    return -mp.tan(mp.pi * s) / (2 * s + 1)


def X(s, lam):
    return mp.gamma(2 * s + lam) * mp.gamma(2 * s + 2) * mp.gamma(s) * mp.gamma(1 - s)


SHAPE = {}


def shape(s):
    key = mp.nstr(s, 25)
    if key not in SHAPE:
        SHAPE[key] = mp.quad(lambda v: mp.e ** (-mp.cosh(v) - s * v),
                             [-16, -8, -4, -2, 0, 2, 4, 8, 16])
    return SHAPE[key]


def tau_grid():
    edges = [0, mp.mpf('0.25'), mp.mpf('0.5'), 1, mp.mpf('1.5'), 2, 3, 4, 6, 8, 12]
    taus, ws = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        m = 4
        h = (b - a) / m
        for j in range(m + 1):
            taus.append(a + j * h)
            ws.append(h / 3 * (1 if j in (0, m) else (4 if j % 2 == 1 else 2)))
    return [-t for t in reversed(taus)] + taus, [w for w in reversed(ws)] + ws


TAUS, WTS = tau_grid()


def G0hat(s, lam, z0, nmax=400):
    tol = mp.mpf(10) ** (-(DPS - 10))
    total = mp.mpf(0)
    Pi = 1 / (lam + 2 * s)
    maxterm = mp.mpf(0)
    n_used = 0
    for n in range(nmax):
        term = Pi * shape(s + mp.mpf(n) / 2) * mp.e ** (-(s + mp.mpf(n) / 2) * z0)
        total += term
        a = abs(term)
        if a > maxterm:
            maxterm = a
        n_used = n
        if n > 6 and a < tol * (abs(total) + a):
            break
        Pi = Pi * c_sym(s + mp.mpf(n + 1) / 2) / (lam + 2 * s + (n + 1))
    return total, maxterm, n_used


FIX = {}


def defect_fix(lam, z0):
    """G0hat(1/2)/X(1/2) by symmetric limit — s=1/2 is a removable zero-x-pole point
    (odd-n ladder terms vanish, even-n have finite c(1)c(3/2)-type limits)."""
    key = (mp.nstr(lam, 20), mp.nstr(z0, 10))
    if key not in FIX:
        d = mp.mpf('1e-8')
        gp, _, _ = G0hat(mp.mpf('0.5') + d, lam, z0)
        gm, _, _ = G0hat(mp.mpf('0.5') - d, lam, z0)
        FIX[key] = ((gp + gm) / 2) / (2 * mp.pi * mp.gamma(1 + lam))
    return FIX[key]


def Ghat(s, lam, z0):
    g0, mt, nu = G0hat(s, lam, z0)
    return g0 - defect_fix(lam, z0) * X(s, lam), mt, nu


def norm_F(z0):
    val = mp.quad(lambda v: mp.e ** (-2 * mp.cosh(v) - 2 * SIGMA * v), [-16, -4, 0, 4, 16])
    return mp.sqrt(val) * mp.e ** (-SIGMA * z0)


def gains_for(lam, z0):
    nG2 = mp.mpf(0)
    nA2 = mp.mpf(0)
    canc = 0.0
    nmx = 0
    for t, w in zip(TAUS, WTS):
        s = SIGMA + 1j * t
        g, mt, nu = Ghat(s, lam, z0)
        fh = shape(s) * mp.e ** (-s * z0)
        nG2 += w * abs(g) ** 2
        nA2 += w * abs(fh / (lam + 2 * s)) ** 2
        if abs(g) > 0 and mt > abs(g):
            canc = max(canc, float(mp.log10(mt / abs(g))))
        nmx = max(nmx, nu)
    nF = norm_F(z0)
    return mp.sqrt(nG2 / (2 * mp.pi)) / nF, mp.sqrt(nA2 / (2 * mp.pi)) / nF, canc, nmx


def line_selfcheck(lam, z0):
    s = SIGMA + 1j * mp.mpf('0.7345')
    gv, _, _ = Ghat(s, lam, z0)
    gv2, _, _ = Ghat(s + mp.mpf('0.5'), lam, z0)
    fh = shape(s) * mp.e ** (-s * z0)
    return abs((lam + 2 * s) * gv - c_sym(s + mp.mpf('0.5')) * gv2 - fh) / abs(fh)


def main():
    LAMS = ([mp.mpf('0.3'), mp.mpf('0.05'), mp.mpf('-0.2'), mp.mpf('-0.5'),
             mp.mpc('0.3', '1.5'), mp.mpf('-0.7')] if DEEP
            else [mp.mpf('0.3'), mp.mpf('-0.5')])
    Z0S = ([mp.mpf('-2'), mp.mpf('-6'), mp.mpf('-10'), mp.mpf('-14')] if DEEP
           else [mp.mpf('-2'), mp.mpf('-6'), mp.mpf('-10')])
    print(f"sigma = {float(SIGMA)}  transport line Re lam = {float(-2 * SIGMA)}  "
          f"dps = {DPS}  ({'deep' if DEEP else 'quick'} mode)")
    ok = True
    resid = line_selfcheck(LAMS[0], Z0S[0])
    print(f"  line-eq self-check: rel resid = {mp.nstr(resid, 3)}")
    ok &= resid < mp.mpf(10) ** (-(DPS - 25))
    t0 = time.time()
    for z0 in Z0S:
        for lam in LAMS:
            g, gA, canc, nmx = gains_for(lam, z0)
            finite = mp.isfinite(g) and g > 0
            ok &= bool(finite)
            print(f"  lam={mp.nstr(lam, 3):>10} z0={float(z0):+3.0f}:  gain={mp.nstr(g, 4):>9}  "
                  f"transport={mp.nstr(gA, 4):>9}  dressing={mp.nstr(g / gA, 4):>8}  "
                  f"[canc 10^{canc:.0f}, n<={nmx}, {time.time() - t0:.0f}s]", flush=True)
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
