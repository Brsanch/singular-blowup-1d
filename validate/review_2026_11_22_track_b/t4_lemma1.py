"""Task 4: Lemma 1 bounds (2.1), (2.2), (2.3), (2.4) on families G = N h with a(G)=0.
Closed-form hhat; functions on the complex lines Im z = y obtained by inverse FFT of closed-form
transforms times e^{-xi y} (no post-hoc e^{rho|xi|} amplification of FFT noise)."""
import math, sys
import numpy as np
from scipy.special import loggamma
from lib import *


def logN(m):
    m = np.asarray(m, dtype=np.complex128)
    s = np.where(m.imag >= 0, 1.0, -1.0)
    return -1j * s * np.pi * m - math.log(2.0) + np.log1p(np.exp(2j * s * np.pi * m)) + loggamma(2.0 * m + 2.0) - math.log(math.pi)


def cot(zeta):
    zeta = np.asarray(zeta, dtype=np.complex128)
    s = np.where(zeta.imag >= 0, 1.0, -1.0)
    e = np.exp(2j * s * np.pi * zeta)  # |e| <= 1
    return -1j * s * (1 + e) / (1 - e)


def Qsym(zeta):
    return cot(zeta) / (zeta + 1.0)


class Fam:
    def __init__(self, name, loghhat):
        self.name, self.lh = name, loghhat

    def logGhat(self, m):
        return logN(m) + self.lh(m)


def fam_pb(p, b, tag=""):
    return Fam(f"h=u^({p})e^(-{b:g}u){tag}", lambda m: math.log(2.0) - (2 * m + p) * math.log(b) + loggamma(2 * m + p))


def fam_gauss(z0, s):
    return Fam(f"h=gauss(z0={z0},s={s})", lambda m: math.log(s * math.sqrt(2 * math.pi)) - m * z0 + 0.5 * m * m * s * s)


g = Grid(-90.0, 90.0, 2**14)
x, xi = g.z, g.xi
PH = np.exp(1j * xi * g.zmin)


def line_fn(logFhat_vals, sig, y):
    """F(x+iy) from Fhat(sig+i xi) (given as log)."""
    v = np.exp(np.clip((logFhat_vals - xi * y).real, -745, 700) + 1j * (logFhat_vals - xi * y).imag)
    return np.exp(sig * x) * np.exp(1j * sig * y) * np.fft.ifft(v * PH) / g.dz


def strip_norm(logFhat_vals, rho):
    return math.sqrt(float(np.sum(np.exp(np.clip(2 * rho * np.abs(xi) + 2 * logFhat_vals.real, -745, 700))) * g.dxi / (2 * np.pi)))


def analyse(F, rho, rhop):
    m1, m2 = 0.25 + 1j * xi, 0.625 + 1j * xi
    lG1, lG2 = F.logGhat(m1), F.logGhat(m2)
    n1, n2 = strip_norm(lG1, rho), strip_norm(lG2, rho)
    nX = n1 + n2
    out = dict(nX=nX, n1=n1, n2=n2)
    # --- (2.1): sup|QG| on the boundary lines y = +-rho (max modulus => sup over the strip)
    z1, z2 = -0.25 + 1j * xi, 0.125 + 1j * xi
    lQ1, lQ2 = np.log(Qsym(z1)) + lG1, np.log(Qsym(z2)) + lG2
    sup = 0.0; supw_pos = 0.0; supw_neg = 0.0; mism = 0.0
    for y in (-rho, rho, 0.0):
        qa = line_fn(lQ1, -0.25, y)
        qb = line_fn(lQ2, 0.125, y)
        core = (x > -40) & (x < 40)
        mism = max(mism, float(np.max(np.abs(qa - qb)[core])) / max(float(np.max(np.abs(qa[core]))), 1e-300))
        sup = max(sup, float(np.max(np.abs(qa[core]))))
        supw_pos = max(supw_pos, float(np.max((np.exp(0.25 * x) * np.abs(qa))[core])))
        supw_neg = max(supw_neg, float(np.max((np.exp(-0.125 * x) * np.abs(qb))[core])))
    out.update(sup=sup, CQ=sup / nX, C2pos=supw_pos / n1, C2neg=supw_neg / n2, two_line_mismatch=mism)
    # --- (2.2): ||P(G)||_{X_rho'} (rho - rho') / ||G||_{X_rho}^2
    nP = []
    for sig, lG in ((0.25, lG1), (0.625, lG2)):
        T = []
        for y in (-rhop, rhop):
            Gy = line_fn(lG1, 0.25, y)
            DGy = line_fn(np.log(m1) + lG1, 0.25, y)
            Qy = line_fn(lQ1, -0.25, y)
            Py = -Qy * DGy + 0.5 * Qy * Gy
            # transform on line sig, times e^{rho'|xi|} for the appropriate sign of xi
            T.append(np.exp(-1j * sig * y) * g.hat(Py, sig))
        Tw = np.where(xi >= 0, T[0], T[1])
        nP.append(math.sqrt(float(np.sum(np.abs(Tw) ** 2) * g.dxi / (2 * np.pi))))
    out.update(nP=sum(nP), CP=sum(nP) * (rho - rhop) / nX**2)
    return out


if __name__ == "__main__":
    # rigorous constants from the proof (Cauchy-Schwarz):  C2(zeta0) = ||Qsym(zeta0+i.)||_L2 / sqrt(2 pi)
    xx = np.linspace(-4000, 4000, 800001)
    for z0 in (-0.25, 0.125):
        c2 = math.sqrt(np.trapezoid(np.abs(Qsym(z0 + 1j * xx)) ** 2, xx) / (2 * np.pi))
        Md = float(np.max(np.abs((z0 + 1j * xx) * Qsym(z0 + 1j * xx))))
        print(f"proof constants on Re zeta={z0}: C2 = ||sym||_2/sqrt(2pi) = {c2:.4f};  sup|zeta*sym| = {Md:.4f}")

    fams = []
    for p in (0.5, 1.0, 3.0, 8.0):
        fams.append(fam_pb(p, 1.0))
    for b in (1e-3, 1e-1, 10.0, 1e3):
        fams.append(fam_pb(3.0, b, tag=" [dilated]"))
    for k in (4.0, 16.0, 40.0):
        fams.append(fam_pb(3.0 + 1j * k, 1.0, tag=" [oscillatory, complex]"))
    for z0, s in ((-12.0, 1.0), (10.0, 1.0), (0.0, 0.25), (-12.0, 0.25), (0.0, 4.0)):
        fams.append(fam_gauss(z0, s))
    print("\n(2.1)/(2.2):  CQ = sup_strip|QG| / ||G||_X ;  C2pos = sup e^{x/4}|QG| / ||G||_{rho,1/4} ; C2neg = sup e^{-x/8}|QG| / ||G||_{rho,5/8};"
          "  CP = ||P(G)||_{X_rho'} (rho-rho') / ||G||_X^2")
    for rho, rhop in ((0.0, 0.0), (0.5, 0.25), (1.0, 0.5), (math.pi / 2, math.pi / 4), (math.pi / 2, 0.98 * math.pi / 2)):
        print(f"\n rho={rho:.4f} rho'={rhop:.4f}")
        mx = dict(CQ=0, C2pos=0, C2neg=0, CP=0)
        for F in fams:
            if rho > 0 and "gauss" in F.name and "s=0.25" in F.name and rho > 1.2:
                continue  # e^{rho^2/(2 s^2)} overflow territory (norm finite but ~1e17); skip
            r = analyse(F, rho, rhop)
            for k_ in mx:
                mx[k_] = max(mx[k_], r[k_])
            print(f"   {F.name:42s} ||G||_X={r['nX']:.3e} sup|QG|={r['sup']:.3e}  CQ={r['CQ']:.4f} C2pos={r['C2pos']:.4f} C2neg={r['C2neg']:.4f}"
                  f"  2-line mismatch={r['two_line_mismatch']:.1e}  CP={r['CP']:.4f}")
        print("   MAX:", {k_: round(v, 4) for k_, v in mx.items()})
