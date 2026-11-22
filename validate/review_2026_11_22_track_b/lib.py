"""Reviewer's own toolkit (written from the conventions in the task, not copied from the gates).

z = ln w, u = e^{-z/2}.  Ghat(m) = int e^{-m z} G dz.
N(m) = cos(pi m) Gamma(2m+2)/pi.
Real-space Q:  (QG)(z) = (e^{-z}/pi) int e^{z'/2} G(z') ln|e^{z-z'}-1| dz'   (no Mellin multiplier).
"""
import math
import numpy as np
from scipy.special import loggamma

ZP0 = -0.9189385332046727  # zeta'(0) = -ln(2 pi)/2
ZP2 = -0.030448457058393270  # zeta'(-2)
ZP4 = 0.007983811450268624  # zeta'(-4)


def Nsym(m):
    """N(m) for complex array m, stable for large |Im m|."""
    m = np.asarray(m, dtype=np.complex128)
    xi = m.imag
    s = np.where(xi >= 0, 1.0, -1.0)
    # cos(pi m) = (1/2) e^{-i s pi m} (1 + e^{2 i s pi m}),  |e^{2 i s pi m}| = e^{-2 pi |xi|}
    logcos = -1j * s * np.pi * m - math.log(2.0) + np.log1p(np.exp(2j * s * np.pi * m))
    return np.exp(logcos + loggamma(2.0 * m + 2.0)) / np.pi


class Grid:
    def __init__(self, zmin, zmax, n):
        self.zmin, self.zmax, self.n = zmin, zmax, n
        self.L = zmax - zmin
        self.dz = self.L / n
        self.z = zmin + np.arange(n) * self.dz
        self.xi = 2 * np.pi * np.fft.fftfreq(n, self.dz)
        self.dxi = 2 * np.pi / self.L
        self.u = np.exp(-0.5 * self.z)
        self.w = np.exp(self.z)

    # weighted multiplier on line Re m = sig:  out = e^{sig z} F^{-1}[ M(sig+i xi) F[e^{-sig z} f] ]
    def mult(self, f, sig, M):
        wg = np.exp(-sig * self.z)
        out = np.fft.ifft(M * np.fft.fft(wg * f)) / wg
        return out

    def hat(self, f, sig):
        """continuous transform samples fhat(sig + i xi_k) (phase included)."""
        return self.dz * np.exp(-1j * self.xi * self.zmin) * np.fft.fft(np.exp(-sig * self.z) * f)


def conj_ops(g, sig=0.25, real=True):
    """Conjugated (Mellin) operators on line sig for G-class, sig-1/2 for Q-class."""
    Ng = Nsym(sig + 1j * g.xi)
    Nr = Nsym(sig - 0.5 + 1j * g.xi)
    mg = sig + 1j * g.xi
    mr = sig - 0.5 + 1j * g.xi
    R = (lambda a: a.real) if real else (lambda a: a)
    ops = dict(
        N=lambda h: R(g.mult(h, sig, Ng)),
        Ninv=lambda G: R(g.mult(G, sig, 1.0 / Ng)),
        Nr=lambda h: R(g.mult(h, sig - 0.5, Nr)),
        D=lambda G: R(g.mult(G, sig, mg)),
        Dr=lambda F: R(g.mult(F, sig - 0.5, mr)),
    )
    ops["Q"] = lambda G: -2.0 * ops["Nr"](g.u * ops["Ninv"](G))
    return ops


def Q_real_matrix(z, order=2):
    """Dense real-space Q on a uniform z grid; G assumed negligible outside the window
    (in the sense that e^{z'/2} G ln|.| is negligible).  Punctured trapezoid with
    Euler-Maclaurin log-singularity corrections:
       int g(d) ln|d| dd ~ h sum' g ln|d| + sum_k 2 zeta'(-2k) h^{2k+1} g^{(2k)}(0)/(2k)!
    Returns (A, C0, C2, C4): QG = A@G + C0*G + C2*(F'') + C4*(F'''') handled by Q_real().
    """
    n = len(z)
    h = z[1] - z[0]
    d = z[:, None] - z[None, :]
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        K = np.where(d > 1.0, d + np.log1p(-np.exp(-np.abs(d))),
                     np.where(d < -1.0, np.log1p(-np.exp(-np.abs(d))), np.log(np.abs(np.expm1(d)))))
    K[np.arange(n), np.arange(n)] = 0.0
    A = (np.exp(-z)[:, None] / np.pi) * K * np.exp(0.5 * z)[None, :] * h
    return A


def Q_real(G, z, A=None):
    """(QG)(z_i) by real-space log-kernel quadrature (no Mellin multiplier)."""
    n = len(z)
    h = z[1] - z[0]
    if A is None:
        A = Q_real_matrix(z)
    # e^{-z} F^{(2k)} = (d/dz + 1)^{2k} P with P = e^{-z/2} G  (avoids e^{z/2} noise amplification)
    P = np.exp(-0.5 * z) * G
    xi = 2 * np.pi * np.fft.fftfreq(n, h)
    Ph = np.fft.fft(P)
    cplx = np.iscomplexobj(G)
    terms = []
    for k in (1, 2, 3):
        t = np.fft.ifft((1.0 + 1j * xi) ** (2 * k) * Ph)
        terms.append(t if cplx else t.real)
    ZP6 = -0.005899759143515937
    corr = (h * P * (math.log(h) + 2 * ZP0) + 2 * ZP2 * h**3 * terms[0] / 2.0
            + 2 * ZP4 * h**5 * terms[1] / 24.0 + 2 * ZP6 * h**7 * terms[2] / 720.0)
    return A @ G + corr / np.pi


def a_quad(G, z):
    """a(G) = -(1/pi) int e^{-z/2} G dz, plain trapezoid on the window."""
    return -(1.0 / np.pi) * np.sum(np.exp(-0.5 * z) * G) * (z[1] - z[0])
