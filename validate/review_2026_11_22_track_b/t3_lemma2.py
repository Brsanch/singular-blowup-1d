"""Task 3: Lemma 2 independently.  Strip norms are computed WITHOUT multiplying FFT output by e^{rho|xi|}:
the function is evaluated on the complex lines Im z = -/+ rho, whose x-FFT *is* fhat(sigma+i xi) e^{rho|xi|}
(for xi >0 / <0).  S_g(t)F = N(e^{-ut} f),  f = N^{-1}F taken on the line 1/4.
"""
import math, sys
import numpy as np
from lib import *

g = Grid(-70.0, 90.0, 2**15)
x, xi = g.z, g.xi
S1, S2 = 0.25, 0.625
absN = {s: np.abs(Nsym(s + 1j * xi)) for s in (S1, S2)}


def T_line(fm, fp, sig, rho):
    """fhat(sig+i xi) e^{rho|xi|} from samples on Im z = -rho (fm) and +rho (fp)."""
    Tm = np.exp(1j * sig * rho) * g.hat(fm, sig)
    Tp = np.exp(-1j * sig * rho) * g.hat(fp, sig)
    return np.where(xi >= 0, Tm, Tp)


def nrm(T, sig):
    return math.sqrt(float(np.sum((absN[sig] * np.abs(T)) ** 2) * g.dxi / (2 * np.pi)))


def phihat_w(sig, rho, beta):
    """beta^{2m-1} pi/cos(pi m) e^{rho|xi|} (analytic continuation; valid on both lines)."""
    m = sig + 1j * xi
    s = np.where(xi >= 0, 1.0, -1.0)
    logcos = -1j * s * np.pi * m - math.log(2.0) + np.log1p(np.exp(2j * s * np.pi * m))
    return np.exp((2 * m - 1) * math.log(beta) + math.log(math.pi) - logcos + rho * np.abs(xi))


def run_case(name, freg, c, beta, rho, ts, fhat_def=None):
    """f = freg + c * u/(u^2+beta^2);  a(F) = -2c/pi.   freg: callable on complex z (pole-free part).
    If fhat_def is given (callable m -> Fhat(m), closed form), f is built from Fhat/N instead and
    a(F) = -(1/pi) Fhat(1/2)."""
    zm, zp = x - 1j * rho, x + 1j * rho
    um, up = np.exp(-0.5 * zm), np.exp(-0.5 * zp)
    if fhat_def is None:
        a = -2.0 * c / math.pi
        phim, phip = um / (um**2 + beta**2), up / (up**2 + beta**2)
        frm, frp = freg(zm), freg(zp)
        TF = {}
        for s in (S1, S2):
            TF[s] = T_line(frm, frp, s, rho) + c * phihat_w(s, rho, beta)
        fm, fp = frm + c * phim, frp + c * phip
    else:
        a = (-(1 / math.pi) * fhat_def(np.array([0.5 + 0j]))[0]).real
        TF = {}
        for s in (S1, S2):
            TF[s] = fhat_def(s + 1j * xi) * np.exp(np.clip(rho * np.abs(xi), 0, 300)) / Nsym(s + 1j * xi)
        # f on Im z = -/+rho by inverse transform on the line 1/4
        fh = fhat_def(S1 + 1j * xi) / Nsym(S1 + 1j * xi)
        ph = np.exp(1j * xi * g.zmin)
        fm = np.exp(S1 * x) * np.exp(-1j * S1 * rho) * np.fft.ifft(fh * np.exp(np.clip(xi * rho, -700, 300)) * ph) / g.dz
        fp = np.exp(S1 * x) * np.exp(+1j * S1 * rho) * np.fft.ifft(fh * np.exp(np.clip(-xi * rho, -700, 300)) * ph) / g.dz
        c, beta = -0.5 * math.pi * a, 1.0  # f~ = f + (pi/2) a phi_*  <=> f = freg + c phi_1
        phim, phip = um / (um**2 + 1), up / (up**2 + 1)
        frm, frp = fm - c * phim, fp - c * phip
    nF1, nF2 = nrm(TF[S1], S1), nrm(TF[S2], S2)
    nX = nF1 + nF2
    print(f"\n[{name}] rho={rho:.4f}  ||F||_1/4={nF1:.4e} ||F||_5/8={nF2:.4e}  a(F)={a:+.4e}  |a|/||F||_X={abs(a)/nX:.4f}")
    print("      t      R1=||SF||_1/4/||F||_1/4   ||SF||_5/8/||F||_X   t^.25||SF||_5/8/|a|   (2.6)ratio   reg-part/||F||_X   rank1 t^.25||N(e^-ut phi)||")
    worst = dict(R1=0, R26=0, reg=0)
    for t in ts:
        em, ep = np.exp(-um * t), np.exp(-up * t)
        n1 = nrm(T_line(em * fm, ep * fp, S1, rho), S1)
        n2 = nrm(T_line(em * fm, ep * fp, S2, rho), S2)
        nreg = nrm(T_line(em * frm, ep * frp, S2, rho), S2)
        nphi = nrm(T_line(em * phim, ep * phip, S2, rho), S2)
        r26 = (n1 + n2) / (nX + t**-0.25 * abs(a))
        s_a = f"{t**0.25*n2/abs(a):10.4f}" if abs(a) > 1e-13 else "       n/a"
        print(f"   {t:8.1e}   {n1/nF1:10.5f}             {n2/nX:12.5f}        {s_a}        {r26:8.4f}     {nreg/nX:10.5f}        {t**0.25*nphi:10.5f}")
        worst["R1"] = max(worst["R1"], n1 / nF1)
        worst["R26"] = max(worst["R26"], r26)
        worst["reg"] = max(worst["reg"], nreg / nX)
    return worst


if __name__ == "__main__":
    ts = (1e-4, 1e-3, 1e-2, 1e-1, 0.3, 1.0)
    # --- rank-one split checks ---
    print("== rank-one checks ==")
    for s in (0.0, 0.25):
        T = g.hat(g.u / (g.u**2 + 1), s)
        ex = math.pi / np.cos(np.pi * (s + 1j * xi[np.abs(xi) < 12]))
        print(f" Mellin(phi_*) on Re m={s}: max|num - pi/cos(pi m)| (|xi|<12) = {np.max(np.abs(T[np.abs(xi)<12]-ex)):.2e}")
    T = g.hat(g.u / (g.u**2 + 1) - 1 / g.u, 0.625)
    sel = np.abs(xi) < 9
    ex = math.pi / np.cos(np.pi * (0.625 + 1j * xi[sel]))
    print(f" Mellin(phi_* - 1/u) on Re m=5/8 vs continuation pi/cos(pi m): {np.max(np.abs(T[sel]-ex)):.2e}")
    for X in (2.0, 4.0, 8.0):
        print(f"   |phihat(5/8+i{X})| e^(pi {X}) = {abs(math.pi/np.cos(np.pi*(0.625+1j*X)))*math.exp(math.pi*X):.6f}  (-> 2 pi = {2*math.pi:.6f})")
    e = 1e-6
    print(" Res_{m=1/2} pi/cos(pi m) =", (e * math.pi / math.cos(math.pi * (0.5 + e)) - e * math.pi / math.cos(math.pi * (0.5 - e))) / 2 * 1 + 0.0,
          " via (m-1/2)*phihat:", e * math.pi / math.cos(math.pi * (0.5 + e)))

    allw = []
    for rho in (0.5, 1.0, math.pi / 2):
        U = lambda z: np.exp(-0.5 * z)
        cases = [
            ("A a=0: f=u^3 e^-u", lambda z: U(z) ** 3 * np.exp(-U(z)), 0.0, 1.0, None),
            ("A a=0: f=u^0.3 e^-u", lambda z: U(z) ** 0.3 * np.exp(-U(z)), 0.0, 1.0, None),
            ("A a=0 edge-concentrated: f=u^3 e^-(u/100)  (w~1e-4)", lambda z: (U(z) / 100) ** 3 * np.exp(-U(z) / 100), 0.0, 1.0, None),
            ("A a=0 far-field: f=(30u)^3 e^-30u", lambda z: (30 * U(z)) ** 3 * np.exp(-30 * U(z)), 0.0, 1.0, None),
            ("C a=0 oscillatory: Re-free u^(3+8i) e^-u", lambda z: U(z) ** (3 + 8j) * np.exp(-U(z)), 0.0, 1.0, None),
            ("F a=0 gaussian in z at z0=-12", lambda z: np.exp(-0.5 * (z + 12) ** 2), 0.0, 1.0, None),
            ("B a!=0 pure rank-1: f=phi_*  (F=u^2 e^-u/2)", lambda z: 0 * z, 1.0, 1.0, None),
            ("B a!=0: f=u^3 e^-u + phi_*", lambda z: U(z) ** 3 * np.exp(-U(z)), 1.0, 1.0, None),
            ("B a!=0 edge: f=phi_beta, beta=100", lambda z: 0 * z, 1.0, 100.0, None),
            ("B a!=0 far: f=phi_beta, beta=0.01", lambda z: 0 * z, 1.0, 0.01, None),
            ("B a!=0 large |a|/norm attempt: f=phi_1 - 0.9*u^1 e^-u", lambda z: -0.9 * U(z) * np.exp(-U(z)), 1.0, 1.0, None),
            ("D a!=0 near strip boundary: f=u^3e^-u/(1+(z/(1.02 rho))^2)+phi_*",
             (lambda r: (lambda z: U(z) ** 3 * np.exp(-U(z)) / (1 + (z / (1.02 * r)) ** 2)))(rho), 1.0, 1.0, None),
            ("D a=0 near strip boundary, edge: f=gauss(z+10)/(1+((z+10)/(1.02 rho))^2)",
             (lambda r: (lambda z: np.exp(-0.5 * (z + 10) ** 2) / (1 + ((z + 10) / (1.02 * r)) ** 2)))(rho), 0.0, 1.0, None),
        ]
        for name, fr, c, beta, _ in (cases if 'full' in sys.argv else []):
            w_ = run_case(name, fr, c, beta, rho, ts)
            allw.append((rho, name, w_))
        # hat-defined F = w/(b+w)^2 (the NM gate's datum; a != 0), b = 2 and b = 1e-4 (edge) and 1e4
        for b in (2.0, 1e-4, 1e4):
            def Fh(m, b=b):
                m = np.asarray(m, dtype=np.complex128)
                sg = np.where(m.imag >= 0, 1.0, -1.0)
                inv_sin = -2j * sg * np.exp(1j * sg * np.pi * m) / (1.0 - np.exp(2j * sg * np.pi * m))
                return np.exp(-(m + 1) * math.log(b)) * np.pi * m * inv_sin
            w_ = run_case(f"H a!=0: F=w/(b+w)^2, b={b:g}", None, None, None, rho, ts, fhat_def=Fh)
            allw.append((rho, f"F=w/(b+w)^2 b={b:g}", w_))
    print("\n== worst ratios over t in [1e-4,1] ==")
    for rho, name, w_ in allw:
        print(f" rho={rho:.3f}  R1max={w_['R1']:8.4f}  (2.6)ratio max={w_['R26']:8.4f}  regpart max={w_['reg']:8.4f}   {name}")
