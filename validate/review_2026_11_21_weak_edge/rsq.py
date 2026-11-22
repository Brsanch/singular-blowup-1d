#!/usr/bin/env python3
"""Reviewer's independent REAL-SPACE operators (no Mellin multiplier anywhere).

H[g](y) = (1/pi) pv int_0^inf g(y')/(y-y') dy',  g = y^{-1/2} G,  z = ln y.
Uses pv int_0^inf y'^{-1/2}/(y-y') dy' = 0 (y>0), so
  H(z) = (1/pi) [ int k(z,z') (G(z')-G(z)) dz' ],  k = e^{z'/2}/(e^z - e^{z'}),
integrand smooth; quadrature = odd-offset ("alternating") rule, step 2dz, which
never touches the diagonal.  Outside the grid window G(z') is taken as 0 and the
-G(z) * int_outside k dz' piece is done in closed form.
(QG)(z) = int_{-inf}^z H(z') e^{z'-z} dz'  (Simpson w/ cubic midpoint interp).
"""
import numpy as np


def H_matrix(z):
    n = len(z)
    dz = z[1] - z[0]
    Y = np.exp(z)
    sY = np.exp(0.5 * z)
    i = np.arange(n)
    odd = ((i[:, None] - i[None, :]) % 2) == 1
    with np.errstate(divide="ignore", invalid="ignore"):
        K = sY[None, :] / (Y[:, None] - Y[None, :])
    K[~odd] = 0.0
    K *= 2.0 * dz
    # virtual odd-offset nodes just outside the grid (G=0 there) so the window never ends at z_i
    jL = np.where(i % 2 == 0, -1, -2)
    jR = np.where((n - i) % 2 == 1, n, n + 1)
    zLn = z[0] + jL * dz
    zRn = z[0] + jR * dz
    kL = 2.0 * dz * np.exp(0.5 * zLn) / (Y - np.exp(zLn))
    kR = 2.0 * dz * np.exp(0.5 * zRn) / (Y - np.exp(zRn))
    tL = np.exp(0.5 * (zLn - dz))
    tR = np.exp(0.5 * (zRn + dz))
    left = np.log(np.abs((sY + tL) / (sY - tL))) / sY
    right = -np.log(np.abs((tR + sY) / (tR - sY))) / sY
    Hm = K.copy()
    Hm[i, i] = -K.sum(axis=1) - kL - kR - (left + right)
    return Hm / np.pi


def EA_apply(Hv, z):
    """q(z_i) = int_{-inf}^{z_i} H e^{z'-z_i} dz'; H assumed const below z[0]."""
    dz = z[1] - z[0]
    n = len(z)
    Hmid = np.empty(n - 1)
    Hmid[1:-1] = (-Hv[:-3] + 9 * Hv[1:-2] + 9 * Hv[2:-1] - Hv[3:]) / 16.0
    Hmid[0] = (5 * Hv[0] + 15 * Hv[1] - 5 * Hv[2] + Hv[3]) / 16.0
    Hmid[-1] = (5 * Hv[-1] + 15 * Hv[-2] - 5 * Hv[-3] + Hv[-4]) / 16.0
    e1 = np.exp(-dz)
    eh = np.exp(-0.5 * dz)
    inc = (dz / 6.0) * (e1 * Hv[:-1] + 4 * eh * Hmid + Hv[1:])
    q = np.empty(n)
    q[0] = Hv[0]
    for k in range(1, n):
        q[k] = e1 * q[k - 1] + inc[k - 1]
    return q


def EA_matrix(z):
    n = len(z)
    M = np.empty((n, n))
    I = np.eye(n)
    for c in range(n):
        M[:, c] = EA_apply(I[:, c], z)
    return M


def D6(G, dz):
    """6th-order central d/dz (one-sided garbage in the outer 3 points -> set nan)."""
    d = np.full_like(G, np.nan)
    d[3:-3] = (-G[:-6] + 9 * G[1:-5] - 45 * G[2:-4] + 45 * G[4:-2] - 9 * G[5:-1] + G[6:]) / (60 * dz)
    return d


def a_quad(G, z, lo=None, hi=None):
    """a = -(1/pi) int e^{-z/2} G dz, trapezoid on window, no tail model."""
    m = np.ones(len(z), bool)
    if lo is not None:
        m &= z >= lo
    if hi is not None:
        m &= z <= hi
    dz = z[1] - z[0]
    f = np.exp(-0.5 * z[m]) * G[m]
    return -(1.0 / np.pi) * (np.sum(f) - 0.5 * (f[0] + f[-1])) * dz


if __name__ == "__main__":
    import mpmath as mp
    b = 2.0
    for (zlo, zhi, n) in ((-24, 24, 1201), (-24, 24, 2401), (-30, 30, 1501), (-18, 18, 901)):
        z = np.linspace(zlo, zhi, n)
        w = np.exp(z)
        Hm = H_matrix(z)
        for name, G, qex in (
            ("w/(b+w)^2", w / (b + w) ** 2,
             np.where(w < 1e-8, -1 / (2 * b**1.5) + 3 * w / (4 * b**2.5),
                      np.log1p(w / b) / (2 * np.sqrt(b) * w) - 1 / (np.sqrt(b) * (b + w)))),
        ):
            Hv = Hm @ G
            q = EA_apply(Hv, z)
            mk = (z > zlo + 4) & (z < zhi - 6)
            print(f"[{zlo},{zhi}] n={n} {name}: max|Q-Qexact| bulk={np.max(np.abs(q-qex)[mk]):.2e}, "
                  f"edge q[0]={q[0]:+.9f} exact a={-1/(2*b**1.5):+.9f}, a_quad={a_quad(G,z):+.9f}")
    # mpmath spot check of the closed form itself (pv integral, no formula reuse)
    mp.mp.dps = 25
    def Hg(y):
        f = lambda t: t**-0.5 * t / (b + t) ** 2
        return (mp.quad(lambda t: (f(t) - f(y)) / (y - t), [0, y / 2, y, 2 * y, 50, mp.inf])
                + f(y) * mp.log(1)) / mp.pi  # pv int_0^inf dt/(y-t) = 0 in the symmetric-limit sense? handled below
    def Qmp(wv):
        # H via substitution t=x^2 : g dt = 2 G(x^2) dx ; pv over x
        Gf = lambda x: x * x / (b + x * x) ** 2
        def Hy(y):
            sy = mp.sqrt(y)
            # (1/pi) pv int_0^inf 2 G(x^2)/(y-x^2) dx ; subtract G(y): pv int_0^inf 2/(y-x^2) dx = 0
            return (1 / mp.pi) * mp.quad(lambda x: 2 * (Gf(x) - Gf(sy)) / (y - x * x), [0, sy / 2, sy, 2 * sy, 20, mp.inf])
        return mp.quad(Hy, [0, wv / 2, wv]) / wv
    for wv in (0.1, 1.0, 7.0):
        qe = mp.log1p(wv / b) / (2 * mp.sqrt(b) * wv) - 1 / (mp.sqrt(b) * (b + wv))
        print(f"mpmath pv check w={wv}: Q={float(Qmp(mp.mpf(wv))):+.12f} closed={float(qe):+.12f}")
