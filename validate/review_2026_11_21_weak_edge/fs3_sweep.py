#!/usr/bin/env python3
"""Task 5: FS-3 (non-analytic bump, eps=0.03) under filter / domain / dt / Nz changes.
Norms are of v = e^{-z/4} h, H^k = sum (1+xi^2)^k |vhat|^2 (same formula as the gate, extended to k=2,3)."""
import sys, math, numpy as np, wk, rsq
np.seterr(all="ignore")
def run(ff, zlo, zhi, Nz, Ns, seed="bump", eps=0.03, spec=False, hist=False):
    ops = wk.Ops(Nz=Nz, zmin=zlo, zmax=zhi, filt_frac=ff)
    z = ops.z
    h0 = wk.seed_bump(ops) if seed == "bump" else wk.seed_u3(ops)
    H = {}
    def rec(j, s, h): H[round(s, 2)] = h.copy()
    wk.etdrk2(ops, h0, eps=eps, Ns=Ns, rec_every=Ns // 12, rec=rec)
    tag = f"filt={'none' if ff is None else round(ff,3)} z=[{zlo},{zhi}] Nz={Nz} dz={ops.dz:.4f} Ns={Ns}"
    hf = H[1.2]
    if not np.all(np.isfinite(hf)):
        bad = min(s for s in H if not np.all(np.isfinite(H[s])))
        print(f"{tag}: NaN from s={bad}"); return
    G = ops.Ng(hf)
    k14 = np.argmin(np.abs(z + 14)); c1 = G[k14] / ops.w[k14]
    aq = rsq.a_quad(G, z, -30, min(20, zhi - 2)) - (1 / math.pi) * c1 * 2 * math.exp(-15)
    ie = np.argmin(np.abs(z + 15))
    n = ops.norms(hf)
    print(f"{tag}: s=1.2 L2={n[0]:.6e} H1={n[1]:.6e} H2={n[2]:.6e} H3={n[3]:.6e} | gate a(z=-15)={ops.Q(G)[ie]:+.2e} a_quad(tail-corr)={aq:+.2e} G/w(-14)={c1:+.4e} 3a(P)={3*rsq.a_quad(ops.P(G), z, -20, min(20, zhi-2)):+.4e}")
    if hist:
        for s in sorted(H):
            n = ops.norms(H[s]); print(f"    s={s:.1f}  L2={n[0]:.5e} H1={n[1]:.5e} H2={n[2]:.5e} H3={n[3]:.5e}")
    if spec:
        for s in (0.0, 0.6, 1.2):
            sp = ops.spectrum(H[s]); ax = np.abs(ops.xi); xm = ops.ximax
            cut = (ff if ff else 1.0) * xm
            bands = [(0, 2), (2, 5), (5, 10), (10, 20), (20, 40), (40, 80), (80, 0.9 * cut), (0.9 * cut, cut), (cut, xm + 1)]
            print(f"    spectrum s={s}: max|vhat| per |xi| band (xi_max={xm:.1f}, cut={cut:.1f}): " +
                  " ".join(f"[{a:.0f},{b:.0f}):{(np.max(sp[(ax>=a)&(ax<b)]) if np.any((ax>=a)&(ax<b)) else 0):.1e}" for a, b in bands))
            p = sp**2
            print(f"      H3-energy share by band: " + " ".join(f"{np.sum(((1+ax**2)**3*p)[(ax>=a)&(ax<b)])/np.sum((1+ax**2)**3*p):.2e}" for a, b in bands))
which = sys.argv[1]
dzb = 72 / 4096
N = lambda lo, hi, f=1.0: 2 * int(round((hi - lo) / (dzb * f) / 2))
if which == "base":
    run(2/3, -36, 36, 4096, 240, spec=True, hist=True)
    for Nz in (1024, 2048, 8192): run(2/3, -36, 36, Nz, 240)
if which == "filt":
    for ff in (None, 0.5, 0.6, 2/3, 0.8): run(ff, -36, 36, 4096, 240, spec=(ff in (None, 0.8)))
if which == "dom":
    for (lo, hi) in ((-48, 48), (-60, 60), (-72, 72), (-96, 96), (-104, 24)):
        for ff in (2/3, None): run(ff, lo, hi, N(lo, hi), 240, spec=(lo == -72 and ff is not None), hist=(lo == -72 and ff is not None))
    for f in (4, 2, 0.5): run(2/3, -72, 72, N(-72, 72, f), 240)
if which == "dt":
    for Ns in (60, 120, 240, 480, 960):
        run(2/3, -36, 36, 4096, Ns); run(2/3, -72, 72, 8192, Ns)
