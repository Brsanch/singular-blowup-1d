#!/usr/bin/env python3
"""(i) where does the filtered scheme blow up on a long-left domain?  (ii) Picard per-iteration diffs
(iii) bump edge slope on a clean domain (iv) a for mutated nonlinearities."""
import sys, math, numpy as np, wk, rsq
np.seterr(all="ignore")
which = sys.argv[1]
if which == "blow":
    for (seed, eps, ff, lo, hi, Nz) in (("bump", 0.03, 2/3, -104, 24, 7282), ("bump", 0.003, 2/3, -104, 24, 7282), ("bump", 3e-5, 2/3, -104, 24, 7282),
                                        ("bump", 0.03, 2/3, -72, 72, 4096), ("bump", 0.003, 2/3, -72, 72, 4096), ("bump", 0.03, None, -72, 72, 4096), ("bump", 0.03, None, -72, 72, 2048)):
        ops = wk.Ops(Nz=Nz, zmin=lo, zmax=hi, filt_frac=ff); z = ops.z
        h0 = wk.seed_bump(ops)
        print(f"-- {seed} eps={eps} filt={ff} z=[{lo},{hi}] Nz={Nz}")
        state = {"n": 0}
        def rec(j, s, h):
            if state["n"] > 14: return
            if j % 10 and np.all(np.isfinite(h)) and np.max(np.abs(ops.w_g*h)) < 1e3*eps: return
            G = ops.Ng(h); P = ops.P(G); v = ops.w_g * P
            if not np.all(np.isfinite(v)):
                print(f"   s={s:.3f}: non-finite"); state["n"] = 99; return
            k = int(np.argmax(np.abs(v))); left = z < -30
            kl = int(np.argmax(np.abs(v) * left))
            state["n"] += 1
            print(f"   s={s:.3f}: max|w_g P|={abs(v[k]):.2e} at z={z[k]:+.1f}; on z<-30: {abs(v[kl]):.2e} at z={z[kl]:+.1f}; "
                  f"max|e^(z/4) QG| on z<-30={np.max(np.abs(ops.w_r[left]**-1*0+ (np.exp(0.25*z)*ops.Q(G))[left])):.2e}; max|e^(-z/4) DG| on z<-30={np.max(np.abs((ops.w_g*ops.Dg(G))[left])):.2e}; ||v||H3={ops.norms(h)[3]:.3e}")
        wk.etdrk2(ops, h0, eps=eps, Ns=240, rec_every=1, rec=rec)
if which == "picard":
    for ff, lo, hi, Nz in ((0.6, -40, 40, 4096), (0.6, -104, 24, 6554), (None, -104, 24, 6554), (0.6, -80, 80, 8192)):
        ops = wk.Ops(Nz=Nz, zmin=lo, zmax=hi, filt_frac=ff); z, u = ops.z, ops.u
        h0 = wk.seed_u3(ops); eps = 0.02; Ns = 240; ds = 1.2 / Ns
        bulk = (z > -12) & (z < 12)
        ex, p1, p2 = wk.etd_weights(u, ds); wt1 = ds * p2; wt0 = ds * (p1 - p2)
        s_grid = np.linspace(0, 1.2, Ns + 1)
        G = np.array([ops.Ng(eps * np.exp(-np.minimum(u * sv, 700)) * h0) for sv in s_grid])
        diffs = []
        for it in range(1, 25):
            fP = np.array([ops.NgI(ops.P(G[j])) for j in range(Ns + 1)])
            Gn = np.zeros_like(G); hc = eps * h0.copy(); Gn[0] = ops.Ng(hc)
            for j in range(Ns):
                hc = ex * hc + wt0 * fP[j] + wt1 * fP[j + 1]; Gn[j + 1] = ops.Ng(hc)
            d = float(np.max(np.abs((Gn - G)[:, bulk]))); diffs.append(d); G = Gn
            if not np.isfinite(d) or d < 1e-12: break
        print(f"Picard filt={ff} z=[{lo},{hi}] Nz={Nz}: diffs = " + " ".join(f"{d:.1e}" for d in diffs))
if which == "mut":
    for (ff, lo, hi, Nz) in ((0.6, -40, 40, 4096), (None, -80, 80, 8192)):
        ops = wk.Ops(Nz=Nz, zmin=lo, zmax=hi, filt_frac=ff); z, w = ops.z, ops.w
        h0 = wk.seed_u3(ops)
        nls = {"true P": ops.P,
               "P with c_adv=-1 (wrong sign)": lambda G: ops.P(G, c_adv=-1.0),
               "P with c_quad=3": lambda G: ops.P(G, c_quad=3.0),
               "G^2 (no Q at all)": lambda G: G * G,
               "0.3*exp(-(z+3)^2/2)*G (linear local)": lambda G: 0.3 * np.exp(-0.5 * (z + 3) ** 2) * G,
               "zero (linear flow)": lambda G: 0.0 * G}
        for name, nl in nls.items():
            H = {}
            def rec(j, s, h): H[j] = h.copy()
            wk.etdrk2(ops, h0, eps=0.02, Ns=240, nl=nl, rec_every=20, rec=rec)
            ie = np.argmin(np.abs(z + 15)); k14 = np.argmin(np.abs(z + 14))
            gate_a = max(abs(ops.Q(ops.Ng(H[j]))[ie]) for j in H if j > 0)
            G = ops.Ng(H[120]); c1 = G[k14] / w[k14]
            aq = rsq.a_quad(G, z, -30, min(20, hi - 2)) - (1 / math.pi) * c1 * 2 * math.exp(-15)
            scale = np.sum(np.exp(-0.5 * z[(z > -30) & (z < 20)]) * np.abs(G[(z > -30) & (z < 20)])) * ops.dz / math.pi
            print(f"filt={ff} z=[{lo},{hi}] nl={name:40s}: gate max|a|(s=0.1..1.2, z=-15)={gate_a:.2e} (thr 5e-6) | a_quad tail-corr (s=0.6)={aq:+.2e} (scale {scale:.1e}) | G/w(-14)={c1:+.4e} vs 3a(nl)={3*rsq.a_quad(nl(G), z, -20, 20):+.4e}")
if which == "slope":
    ops = wk.Ops(Nz=10922, zmin=-96, zmax=96, filt_frac=None); z, w = ops.z, ops.w
    h = wk.etdrk2(ops, wk.seed_bump(ops), eps=0.03, Ns=240)
    G = ops.Ng(h); t = 3 * rsq.a_quad(ops.P(G), z, -20, 20)
    print("bump clean domain, s=1.2: 3a(P)=%.5e ;  G/w at z=" % t + " ".join(f"{zz}:{G[np.argmin(np.abs(z-zz))]/w[np.argmin(np.abs(z-zz))]:+.4e}" for zz in (-8, -9, -10, -11, -12, -13, -14, -16)))
    print("   u^2 h at z=" + " ".join(f"{zz}:{(ops.u**2*h)[np.argmin(np.abs(z-zz))]:+.5e}" for zz in (-10, -14, -20, -30)) + f"  target -pi/2 a(P)={-math.pi/6*t:+.5e}")
