#!/usr/bin/env python3
"""Reviewer's parametrised copies of the gates' conjugated solvers (same formulas as
validate/edge_pure_step_nonlinear_manifold.py and ..._finite_sobolev.py), with the
filter fraction, domain, Nz, Ns and the nonlinearity exposed."""
import math
import numpy as np
from scipy.special import loggamma


def build_Nm(xi, sig):
    mm = sig + 1j * xi
    sgn = np.where(xi >= 0, 1.0, -1.0)
    lc = (np.pi * np.abs(xi) - np.log(2.0) - 1j * np.pi * sig * sgn
          + np.log1p(np.exp(-2.0 * np.pi * np.abs(xi) + 2j * np.pi * sig * sgn)))
    return np.exp(lc + loggamma(2.0 * mm + 2.0) - np.log(np.pi))


class Ops:
    def __init__(self, Nz=4096, zmin=-40.0, zmax=40.0, sigma=0.25, filt_frac=0.6):
        self.Nz, self.zmin, self.zmax, self.sigma = Nz, zmin, zmax, sigma
        self.L = zmax - zmin
        self.dz = self.L / Nz
        self.z = zmin + np.arange(Nz) * self.dz
        self.w = np.exp(self.z)
        self.u = np.exp(-0.5 * self.z)
        self.xi = 2.0 * np.pi * np.fft.fftfreq(Nz, self.dz)
        self.Nm_g = build_Nm(self.xi, sigma)
        self.Nm_r = build_Nm(self.xi, sigma - 0.5)
        self.w_g = np.exp(-sigma * self.z)
        self.w_r = np.exp(-(sigma - 0.5) * self.z)
        self.m_g = sigma + 1j * self.xi
        self.ximax = np.max(np.abs(self.xi))
        if filt_frac is None:
            self.filt = np.ones(Nz)
        else:
            self.filt = np.where(np.abs(self.xi) <= filt_frac * self.ximax, 1.0, 0.0)

    def Ng(self, h):
        return np.fft.ifft((self.Nm_g * self.filt) * np.fft.fft(self.w_g * h)).real / self.w_g

    def NgI(self, G):
        return np.fft.ifft((self.filt / self.Nm_g) * np.fft.fft(self.w_g * G)).real / self.w_g

    def Nr(self, h):
        return np.fft.ifft((self.Nm_r * self.filt) * np.fft.fft(self.w_r * h)).real / self.w_r

    def Q(self, G):
        return -2.0 * self.Nr(self.u * self.NgI(G))

    def Dg(self, G):
        return np.fft.ifft((self.m_g * self.filt) * np.fft.fft(self.w_g * G)).real / self.w_g

    def P(self, G, c_adv=1.0, c_quad=0.5):
        QG = self.Q(G)
        return -c_adv * QG * self.Dg(G) + c_quad * QG * G

    def norms(self, h, kmax=3):
        vhat = np.fft.fft(self.w_g * h) / self.Nz
        p = np.abs(vhat) ** 2
        return [math.sqrt(float(np.sum((1.0 + self.xi**2) ** k * p)) * self.L) for k in range(kmax + 1)]

    def spectrum(self, h):
        vhat = np.fft.fft(self.w_g * h) / self.Nz
        return np.abs(vhat)


def etd_weights(u, ds):
    x = u * ds
    ex = np.exp(-x)
    p1 = np.where(x < 1e-4, 1.0 - 0.5 * x + x**2 / 6.0 - x**3 / 24.0, (1.0 - ex) / np.where(x < 1e-4, 1.0, x))
    p2 = np.where(x < 1e-4, 0.5 - x / 6.0 + x**2 / 24.0 - x**3 / 120.0, (1.0 - p1) / np.where(x < 1e-4, 1.0, x))
    return ex, p1, p2


def seed_u3(ops):
    u, z = ops.u, ops.z
    h0 = np.where(u < 200.0, (u**3) * np.exp(-np.minimum(u, 700.0)), 0.0)
    G0 = ops.Ng(h0)
    bulk = (z > -12.0) & (z < 12.0)
    return h0 / float(np.max(np.abs(G0[bulk])))


def seed_bump(ops):
    u, z = ops.u, ops.z
    inside = (u > 0.5) & (u < 4.0)
    t_u = np.where(inside, (u - 0.5) / 3.5, 0.5)
    h0 = np.where(inside, np.exp(-1.0 / (t_u * (1.0 - t_u))), 0.0)
    G0 = ops.Ng(h0)
    bulk = (z > -12.0) & (z < 12.0)
    return h0 / float(np.max(np.abs(G0[bulk])))


def picard(ops, h0, eps=0.02, S_max=1.2, Ns=240, nl=None, maxit=25, tol=1e-12):
    """Gate PM-2 algorithm verbatim (trapezoid-ETD Picard). Returns s_grid, G[j], h[j], iters, diff."""
    if nl is None:
        nl = ops.P
    z, u = ops.z, ops.u
    bulk = (z > -12.0) & (z < 12.0)
    s_grid = np.linspace(0.0, S_max, Ns + 1)
    ds = S_max / Ns
    ex, p1, p2 = etd_weights(u, ds)
    wt1 = ds * p2
    wt0 = ds * (p1 - p2)
    G = np.array([ops.Ng(eps * np.exp(-np.minimum(u * sv, 700.0)) * h0) for sv in s_grid])
    hh = None
    for it in range(1, maxit):
        fP = np.array([ops.NgI(nl(G[j])) for j in range(Ns + 1)])
        G_new = np.zeros_like(G)
        hh = np.zeros_like(G)
        hc = eps * h0.copy()
        hh[0] = hc
        G_new[0] = ops.Ng(hc)
        for j in range(Ns):
            hc = ex * hc + wt0 * fP[j] + wt1 * fP[j + 1]
            hh[j + 1] = hc
            G_new[j + 1] = ops.Ng(hc)
        diff = float(np.max(np.abs((G_new - G)[:, bulk])))
        G = G_new
        if diff < tol:
            break
    return s_grid, G, hh, it, diff


def etdrk2(ops, h0, eps=0.03, S_max=1.2, Ns=240, nl=None, rec_every=20, rec=None):
    """Gate FS-3 algorithm verbatim (ETD-RK2). rec(j, s, h) callback."""
    if nl is None:
        nl = ops.P
    u = ops.u
    ds = S_max / Ns
    ex, p1, p2 = etd_weights(u, ds)
    f = lambda hc: ops.NgI(nl(ops.Ng(hc)))
    h = eps * h0.copy()
    if rec:
        rec(0, 0.0, h)
    for j in range(Ns):
        f0 = f(h)
        hp = ex * h + (ds * p1) * f0
        f1 = f(hp)
        h = ex * h + ds * (p1 - p2) * f0 + ds * p2 * f1
        if rec and ((j + 1) % rec_every == 0):
            rec(j + 1, (j + 1) * ds, h)
    return h
