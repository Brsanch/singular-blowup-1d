"""Task 4 (cont.): identity (2.3) both ways with a REAL-SPACE Q; Mellin symbol of Q vs real-space Q;
zero-derivative-loss bound (2.4) at rho = 0 on oscillatory families."""
import math
import numpy as np
from lib import *
from t4_lemma1 import Qsym


def norms0(g, G):
    out = []
    for s in (0.25, 0.625):
        out.append(math.sqrt(float(np.sum(np.abs(g.hat(G, s)) ** 2) * g.dxi / (2 * np.pi))))
    return out


def study(g, h, label, use_real=True):
    ops = conj_ops(g)
    G = ops["N"](h)
    z = g.z
    Qc = ops["Q"](G)
    core = (z > -30) & (z < 30)
    msg = ""
    if use_real:
        Qr = Q_real(G, z)
        # Mellin-symbol Q on the two lines
        Qm1 = g.mult(G, 0.25, 1.0).real  # placeholder to keep shapes
        Gh1 = np.fft.fft(np.exp(-0.25 * z) * G)
        Qm1 = (np.fft.ifft(Qsym(-0.25 + 1j * g.xi) * Gh1) * np.exp(-0.25 * z)).real
        Gh2 = np.fft.fft(np.exp(-0.625 * z) * G)
        Qm2 = (np.fft.ifft(Qsym(0.125 + 1j * g.xi) * Gh2) * np.exp(0.125 * z)).real
        c2 = (z > -15) & (z < 30)
        msg = (f" | max|Q_real - Q_conj|={np.max(np.abs(Qr-Qc)[core]):.1e}  |Q_real - Q_sym(line -1/4)|={np.max(np.abs(Qr-Qm1)[core]):.1e}"
               f"  |Q_real - Q_sym(line 1/8)|(z>-15)={np.max(np.abs(Qr-Qm2)[c2]):.1e}  (max|Q|={np.max(np.abs(Qr[core])):.2e})")
        Quse = Qr
    else:
        Quse = Qc
    DG = ops["D"](G)
    DQ = ops["Dr"](Quse)
    P = -Quse * DG + 0.5 * Quse * G
    wgt = np.exp(-0.5 * z)
    mk = (z > -40) & (z < 40)
    a_dir = -(1 / math.pi) * np.sum((wgt * P)[mk]) * g.dz
    a_ibp = -(1 / math.pi) * np.sum((wgt * DQ * G)[mk]) * g.dz
    # second-order finite-difference D (no spectral operator at all) as a third path
    DGfd = np.gradient(G, g.dz, edge_order=2)
    a_fd = -(1 / math.pi) * np.sum((wgt * (-Quse * DGfd + 0.5 * Quse * G))[mk]) * g.dz
    n1, n2 = norms0(g, G)
    d1, d2 = norms0(g, DG)
    aG = a_quad(G, z)
    Ra = math.pi * abs(a_dir) / (n1**2 + n2**2)
    Rloss = abs(a_dir) / ((n1 + n2) * (d1 + d2))
    print(f"{label:34s} a(G)={aG:+.1e} a(P) direct={a_dir:+.10e} ibp={a_ibp:+.10e} (rel {abs(a_dir-a_ibp)/abs(a_dir):.1e}; FD-path rel {abs(a_fd-a_dir)/abs(a_dir):.1e})"
          f"  ||G||_0X={n1+n2:.3e} ||DG||_0X={d1+d2:.3e}  pi|a(P)|/(n1^2+n2^2)={Ra:.4f}  |a(P)|/(||G|| ||DG||)={Rloss:.2e}{msg}")
    return Ra


if __name__ == "__main__":
    g = Grid(-48.0, 48.0, 4096)
    u = g.u
    print("== real-space Q (n=4096) ==")
    study(g, u**3 * np.exp(-u), "h=u^3 e^-u")
    study(g, u**1 * np.exp(-u), "h=u e^-u")
    study(g, u**3 * np.exp(-u) * np.cos(4 * np.log(u)), "h=u^3 e^-u cos(4 ln u)")
    study(g, u**3 * np.exp(-u) * np.cos(12 * np.log(u)), "h=u^3 e^-u cos(12 ln u)")
    study(g, np.exp(-0.5 * (g.z + 8) ** 2), "h=gauss(z+8)")
    study(g, np.exp(-0.5 * (g.z - 8) ** 2), "h=gauss(z-8)")
    study(g, np.exp(-0.5 * (g.z) ** 2 / 0.3**2), "h=gauss(z/0.3)")
    print("== conjugated Q on fine grid (n=16384), oscillation sweep: is a(P(G)) controlled by ||G||^2 with no derivative? ==")
    g2 = Grid(-48.0, 48.0, 16384)
    u = g2.u
    mx = 0
    for k in (0, 2, 4, 8, 16, 32, 64, 128):
        mx = max(mx, study(g2, u**3 * np.exp(-u) * np.cos(k * np.log(u)), f"h=u^3 e^-u cos({k} ln u)", use_real=False))
    for s in (2.0, 1.0, 0.5, 0.25, 0.1, 0.05):
        mx = max(mx, study(g2, np.exp(-0.5 * (g2.z) ** 2 / s**2), f"h=gauss(z/{s})", use_real=False))
    for s in (0.5, 0.1):
        mx = max(mx, study(g2, np.exp(-0.5 * (g2.z + 10) ** 2 / s**2), f"h=gauss((z+10)/{s})", use_real=False))
        mx = max(mx, study(g2, np.exp(-0.5 * (g2.z - 10) ** 2 / s**2), f"h=gauss((z-10)/{s})", use_real=False))
    print("max pi|a(P)|/(||G||_{1/4}^2+||G||_{5/8}^2) =", mx, " (proof bound: <= 1)")
