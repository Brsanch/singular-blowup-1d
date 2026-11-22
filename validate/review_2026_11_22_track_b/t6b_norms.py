"""Task 6 (b): the theorem's solution. Track both line norms over s in [0,1.5] at several analyticity
radii, measured analyticity radius, grid-convergence, c1 = 3 a(P(G)) (real-space Q), a(G), horizon ODE."""
import math, sys
import numpy as np
from lib import *
from solver import *
from scipy.interpolate import CubicSpline

RHOS = (0.0, 0.25, 0.5, math.pi / 4, 1.0, math.pi / 2)
P_, B_ = 3.0, 1.0
sc = G0_maxnorm_scale(P_, B_)
S = 1.5
NSNAP = 15


def smstep(t):
    t = np.clip(t, 0.0, 1.0)
    a = np.where(t > 0, np.exp(-1.0 / np.maximum(t, 1e-300)), 0.0)
    b = np.where(t < 1, np.exp(-1.0 / np.maximum(1 - t, 1e-300)), 0.0)
    return a / (a + b)


def spec14(g, h):
    """|Ghat(1/4 + i xi)|, xi >= 0, from h on the real line (Ghat = N hhat; line 1/4 is noise-safe)."""
    return (np.abs(Nsym(0.25 + 1j * g.xi)) * np.abs(g.hat(h, 0.25)))[: g.n // 2]


def spec58(g, h, ext=4, zlo=-22.0, zhi=-14.0, ret_c=False):
    """|Ghat(5/8 + i xi)|, xi >= 0.  e^{-5z/8} h has the slow edge tail c e^{3z/8} (h ~ c2/u^2 + c3/u^3 + c4/u^4),
    which is both noise-amplified and truncated on the solver grid.  Replace h below z ~ -18 by the fitted
    tail (smooth blend on [zlo,zhi]) on a grid extended to the left by (ext-1) L, then FFT there."""
    z = g.z
    u = np.exp(-0.5 * z)
    m = (z > zlo) & (z < zhi)
    t = 1.0 / u[m]
    M = np.vstack([np.ones(m.sum()), t, t * t]).T
    co = np.linalg.lstsq(M, (h * u * u)[m], rcond=None)[0]
    ne = ext * g.n
    ze = g.zmax - (np.arange(ne)[::-1] + 1) * g.dz
    te = np.exp(0.5 * ze)  # 1/u
    tail = te**2 * (co[0] + co[1] * te + co[2] * te**2)
    chi = smstep((ze - zlo) / (zhi - zlo))
    he = np.zeros(ne)
    he[-g.n:] = h
    hd = chi * he + (1 - chi) * tail
    F = np.fft.fft(np.exp(-0.625 * ze) * hd) * g.dz
    A = np.abs(F[::ext])[: g.n // 2]
    out = np.abs(Nsym(0.625 + 1j * g.xi[: g.n // 2])) * A
    return (out, co) if ret_c else out


def wnorm(g, A, rho, floor=1e-10):
    """sqrt((1/2pi) int e^{2 rho|xi|}|Ghat|^2) from |Ghat| on xi>=0 (real G), cut where |Ghat| < floor*max.
    Returns (norm, resolved?) ; resolved = weighted integrand has decayed by 1e-6 at the cutoff."""
    xi = g.xi[: g.n // 2]
    mx = A.max(); kp = int(np.argmax(A))
    below = np.where((A < floor * mx) & (np.arange(len(A)) > kp))[0]
    kc = below[0] if len(below) else len(A)
    integ = np.exp(2 * rho * xi[:kc]) * A[:kc] ** 2
    val = 2 * np.sum(integ) - integ[0]
    return math.sqrt(float(val) * g.dxi / (2 * np.pi)), bool(integ[kc - 1] < 1e-6 * integ.max())


def app_rate(g, A):
    """apparent exponential decay rate of |Ghat| between relative levels 1e-6 and 1e-10 (past the peak)."""
    xi = g.xi[: g.n // 2]
    mx = A.max(); kp = int(np.argmax(A))
    idx = np.arange(len(A))
    k1 = np.where((A < 1e-6 * mx) & (idx > kp))[0]
    k2 = np.where((A < 1e-10 * mx) & (idx > kp))[0]
    if not len(k1) or not len(k2) or k2[0] <= k1[0]:
        return float("nan")
    return float(math.log(A[k1[0]] / A[k2[0]]) / (xi[k2[0]] - xi[k1[0]]))


def Xnorm_h(g, h, rho):
    a1, r1 = wnorm(g, spec14(g, h), rho)
    a2, r2 = wnorm(g, spec58(g, h), rho)
    return a1 + a2, (r1 and r2), a1, a2


def VG(g, G, r):
    z = g.z
    arg = r * r * np.exp(-z)
    ln_term = np.where(arg > 1e6, 2.0 * math.log(r) - z + np.log1p(1.0 / np.maximum(arg, 1e-300)), np.log1p(arg))
    return float(-(1.0 / (math.pi * r)) * np.sum(np.exp(0.5 * z) * ln_term * G) * g.dz)


def main(A):
    print(f"\n================ A = max|G_0| = {A}  (h0 = u^{P_} e^(-{B_}u)) ================")
    thm = {}
    for rho0 in (0.5, 1.0, math.pi / 2):
        nx = sum(strip_norms_closed_form_h0(P_, B_, rho0, sg) for sg in (0.25, 0.625)) / sc
        thm[rho0] = A * nx
        print(f"  theorem epsilon = ||eps G_0||_X_rho0  at rho0={rho0:.4f}: {A*nx:.4e}")
    runs = {}
    for n, nst in ((4096, 600), (8192, 1200)):
        s = Solver(-64.0, 64.0, n, sponge=10.0)
        H, bad = s.run(A * h0_family(s, P_, B_) / sc, S, nst, nst // NSNAP)
        runs[n] = (s, H, bad)
        if bad:
            print(f"  n={n}: NON-FINITE at s={bad:.4f}")
    s, H, bad = runs[8192]
    s4, H4, bad4 = runs[4096]
    g = s.g; z = g.z
    core = (s4.g.z > -12) & (s4.g.z < 12)
    Aq = Q_real_matrix(s4.g.z)
    print("   s    max|G|  gridconv(4096 vs 8192)  rho*_1/4 rho*_5/8 [xi_cut] | ||G(s)||_X_rho/||G(0)||_X_rho for rho=0,.25,.5,pi/4,1,pi/2 (u=unresolved/inf)"
          " | thm ratio ||G(s)||_X(rho0/2)/||G(0)||_X(rho0), rho0=.5,1,pi/2 | a(G) quad | c1 fit | 3a(P) realQ | 3a(P)/c1-1")
    n0 = None
    Gs = []
    for j in range(len(H)):
        G = s.N(H[j]); Gs.append(G)
        if j < len(H4):
            G4 = s4.N(H4[j])
            gc = float(np.max(np.abs(CubicSpline(z, G)(s4.g.z[core]) - G4[core])))
        else:
            gc = float("nan")
        r1 = app_rate(g, spec14(g, H[j])); r2 = app_rate(g, spec58(g, H[j])); xc = 0.0
        nr = [Xnorm_h(g, H[j], rho) for rho in RHOS]
        if n0 is None:
            n0 = [v[0] for v in nr]
            cf = [sum(strip_norms_closed_form_h0(P_, B_, rho, sg) for sg in (0.25, 0.625)) * A / sc for rho in RHOS]
            print("   check s=0 norms vs closed form (rel err):", " ".join(f"{n0[i]/cf[i]-1:+.1e}" for i in range(len(RHOS))))
        ratios = " ".join(f"{v[0]/n0[i]:7.4f}{'' if v[1] else 'u'}" for i, v in enumerate(nr))
        thmr = " ".join(f"{nr[i1][0]/n0[i0]:8.5f}{'' if nr[i1][1] else 'u'}" for i1, i0 in ((1, 2), (2, 4), (3, 5)))
        # edge slope and a(P) with the real-space Q (on the 4096 grid)
        extra = ""
        if j < len(H4):
            zz = s4.g.z
            c1 = float(np.mean((G4 / s4.g.w)[(zz > -18) & (zz < -15)]))
            Qr = Q_real(G4, zz, Aq)
            P = -Qr * s4.D(G4) + 0.5 * Qr * G4
            mk = (zz > -40) & (zz < 40)
            aP = float(-(1 / math.pi) * np.sum((np.exp(-0.5 * zz) * P)[mk]) * s4.g.dz)
            aG = float(-(1 / math.pi) * np.sum((np.exp(-0.5 * zz) * G4)[mk]) * s4.g.dz)
            extra = f" | {aG:+.1e} | {c1:+.6e} | {3*aP:+.6e} | {3*aP/c1-1 if abs(c1)>1e-300 else float('nan'):+.1e}"
        print(f"  {j*S/NSNAP:4.2f} {np.max(np.abs(G[(z>-30)&(z<30)])):.3e}  {gc:.1e}   {r1:6.3f} {r2:6.3f} [{xc:5.1f}] | {ratios} | {thmr}{extra}")
    if bad:
        return
    # horizon ODE dr/ds = -1 + V/2 with RK4 over stored snapshots (cubic-in-s interpolation of G not needed:
    # re-run the solver in lockstep with fine steps)
    s = Solver(-64.0, 64.0, 4096, sponge=10.0)
    nst = 1500
    s.setup(S / nst)
    h = A * h0_family(s, P_, B_) / sc
    for mut, lab in ((lambda v: v, "true V_G"), (lambda v: -v, "MUT -V_G"), (lambda v: 3 * v, "MUT 3V_G")):
        s.setup(S / nst)
        h = A * h0_family(s, P_, B_) / sc
        r = 1.0; sinf = None; Vmax = 0.0
        G = s.N(h)
        for k in range(nst):
            # Heun with G at both ends
            v0 = mut(VG(s.g, G, r)); Vmax = max(Vmax, abs(v0))
            k1 = -1 + 0.5 * v0
            h = s.step(h); Gn = s.N(h)
            rp = r + s.ds * k1
            if rp <= 0:
                sinf = k * s.ds - r / k1; slope = k1; break
            k2 = -1 + 0.5 * mut(VG(s.g, Gn, rp))
            rn = r + 0.5 * s.ds * (k1 + k2)
            if rn <= 0:
                sinf = k * s.ds - r / k1; slope = k1; break
            r, G = rn, Gn
        print(f"  horizon ({lab}): s_inf = {sinf if sinf is None else round(sinf,6)}  dr/ds at crossing = {slope if sinf else None}  max|V_G| along path = {Vmax:.3e}")


if __name__ == "__main__":
    for a in sys.argv[1:]:
        main(float(a))
