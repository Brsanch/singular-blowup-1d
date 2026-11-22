"""Reviewer's own integrator for the s-clock pure-step equation in conjugated form
    h_s + u h = N^{-1} P(N h),   P(G) = -(QG) DG + (1/2)(QG) G,   QG = -2 N_r(u h),
ETDRK4 (Cox-Matthews) with pointwise phi-functions (the linear part is multiplication by -u).
Works on the real line (y=0) or on a complex line Im z = y (all operators are x-multipliers
or pointwise products of holomorphic functions, so the same equation holds on the shifted line).
The real-space Q of lib.py is used only for verification (t6_*.py), never inside the solver.
"""
import math
import numpy as np
from lib import *


def phis(x):
    """phi1, phi2, phi3 of complex/real array x (Re x <= 0)."""
    x = np.asarray(x, dtype=np.complex128)
    small = np.abs(x) < 0.2
    xs = np.where(small, x, 1.0)
    # Taylor: phi_k = sum_{j>=0} x^j/(j+k)!
    def tay(k):
        s = np.zeros_like(xs)
        for j in range(18, -1, -1):
            s = s * xs + 1.0 / math.factorial(j + k)
        return s
    xl = np.where(small, 1.0, x)
    e = np.exp(xl)
    p1 = (e - 1) / xl
    p2 = (e - 1 - xl) / xl**2
    p3 = (e - 1 - xl - xl**2 / 2) / xl**3
    return (np.where(small, tay(1), p1), np.where(small, tay(2), p2), np.where(small, tay(3), p3))


class Solver:
    def __init__(self, zmin, zmax, n, y=0.0, sig=0.25, filt=None, sponge=None):
        self.g = g = Grid(zmin, zmax, n)
        self.y = y
        self.real = (y == 0.0)
        self.zc = g.z + 1j * y if not self.real else g.z
        self.u = np.exp(-0.5 * self.zc)
        self.sig = sig
        xi = g.xi
        self.Ng = Nsym(sig + 1j * xi)
        self.Nr = Nsym(sig - 0.5 + 1j * xi)
        self.mg = sig + 1j * xi
        self.wg = np.exp(-sig * g.z)
        self.wr = np.exp(-(sig - 0.5) * g.z)
        self.F = np.ones(n) if filt is None else (np.abs(xi) <= filt * np.max(np.abs(xi))).astype(float)
        # optional end-sponge on the nonlinear forcing (kills the periodic seam where u jumps by e^{L/2})
        self.mask = None
        if sponge:
            def sm(t):  # smooth step 0->1 on [0,1]
                t = np.clip(t, 0.0, 1.0)
                a = np.where(t > 0, np.exp(-1.0 / np.maximum(t, 1e-300)), 0.0)
                b = np.where(t < 1, np.exp(-1.0 / np.maximum(1 - t, 1e-300)), 0.0)
                return a / (a + b)
            self.mask = sm((g.z - zmin) / sponge) * sm((zmax - g.z) / sponge)

    def _R(self, a):
        return a.real if self.real else a

    def N(self, h):
        return self._R(np.fft.ifft(self.F * self.Ng * np.fft.fft(self.wg * h)) / self.wg)

    def Ninv(self, G):
        return self._R(np.fft.ifft(self.F / self.Ng * np.fft.fft(self.wg * G)) / self.wg)

    def Nrr(self, h):
        return self._R(np.fft.ifft(self.F * self.Nr * np.fft.fft(self.wr * h)) / self.wr)

    def D(self, G):
        return self._R(np.fft.ifft(self.F * self.mg * np.fft.fft(self.wg * G)) / self.wg)

    def Q_of_h(self, h):
        return -2.0 * self.Nrr(self.u * h)

    def nl(self, h):
        G = self.N(h)
        q = self.Q_of_h(h)
        f = self.Ninv(-q * self.D(G) + 0.5 * q * G)
        return f if self.mask is None else f * self.mask

    def setup(self, ds):
        x = -self.u * ds
        self.ds = ds
        self.E = np.exp(x)
        self.E2 = np.exp(x / 2)
        p1h, _, _ = phis(x / 2)
        p1, p2, p3 = phis(x)
        cast = (lambda a: a.real) if self.real else (lambda a: a)
        self.E, self.E2 = cast(self.E), cast(self.E2)
        self.a1 = cast(0.5 * ds * p1h)
        self.f1 = cast(ds * (p1 - 3 * p2 + 4 * p3))
        self.f2 = cast(ds * (2 * p2 - 4 * p3))
        self.f3 = cast(ds * (-p2 + 4 * p3))

    def step(self, h):
        Nh = self.nl(h)
        a = self.E2 * h + self.a1 * Nh
        Na = self.nl(a)
        b = self.E2 * h + self.a1 * Na
        Nb = self.nl(b)
        c = self.E2 * a + self.a1 * (2 * Nb - Nh)
        Nc = self.nl(c)
        return self.E * h + self.f1 * Nh + self.f2 * (Na + Nb) + self.f3 * Nc

    def run(self, h0, S, nsteps, store_every=1):
        self.setup(S / nsteps)
        h = h0.copy()
        out = [h.copy()]
        for k in range(1, nsteps + 1):
            h = self.step(h)
            if not np.all(np.isfinite(h)):
                print(f"   !! non-finite at step {k} (s={k*self.ds:.4f})")
                return np.array(out), k * self.ds
            if k % store_every == 0:
                out.append(h.copy())
        return np.array(out), None


def h0_family(sol, p=3.0, b=1.0):
    return sol.u**p * np.exp(-b * sol.u)


def G0_maxnorm_scale(p=3.0, b=1.0):
    """max |N h0| on a reference real grid."""
    s = Solver(-64.0, 64.0, 4096)
    return float(np.max(np.abs(s.N(h0_family(s, p, b)))))


def strip_norms_closed_form_h0(p, b, rho, sig):
    """|| N(u^p e^{-bu}) ||_{rho,sig} from the closed-form transform 2 b^{-(2m+p)} Gamma(2m+p) N(m)."""
    from scipy.special import loggamma
    xi = np.linspace(-200, 200, 400001)
    m = sig + 1j * xi
    lg = math.log(2.0) - (2 * m + p) * math.log(b) + loggamma(2 * m + p)
    val = np.abs(Nsym(m)) ** 2 * np.exp(np.clip(2 * lg.real + 2 * rho * np.abs(xi), -745, 700))
    return math.sqrt(float(np.trapezoid(val, xi)) / (2 * np.pi))


def weighted_norm_real_line(g, G, sig, rho, floor=1e-12, return_diag=False):
    """||G||_{rho,sig} from real-line samples, with a noise-aware cutoff in xi:
    integrate e^{2 rho|xi|}|Ghat|^2 only up to the first |xi| where |Ghat| < floor * max|Ghat|.
    Also returns the measured exponential decay rate of |Ghat| (analyticity radius estimate)."""
    Gh = np.abs(g.hat(G, sig))
    pos = g.xi >= 0
    xi = g.xi[pos][: g.n // 2]
    A = np.maximum(Gh[pos][: g.n // 2], Gh[::-1][: g.n // 2] * 0)  # real G: symmetric
    mx = A.max()
    below = np.where(A < floor * mx)[0]
    kc = below[0] if len(below) else len(A)
    integ = np.exp(2 * rho * xi[:kc]) * A[:kc] ** 2
    val = 2 * np.sum(integ) - integ[0]
    nrm = math.sqrt(float(val) * g.dxi / (2 * np.pi))
    # decay-rate fit of ln|Ghat| between 1e-3*max and 1e-10*max (past the peak)
    kp = int(np.argmax(A))
    sel = np.arange(kp, kc)
    sel = sel[(A[sel] < 1e-3 * mx) & (A[sel] > 1e-10 * mx)]
    rate = float("nan")
    if len(sel) > 5:
        rate = -np.polyfit(xi[sel], np.log(A[sel]), 1)[0]
    # is the weighted integrand still decaying at the cutoff?
    resolved = bool(kc > 3 and integ[kc - 1] < 1e-6 * integ.max())
    if return_diag:
        return nrm, rate, resolved, xi[kc - 1]
    return nrm, rate, resolved
