#!/usr/bin/env python3
"""(A) V_G-mutation of the gates' physical-time checks (PM-3, NM-3c): do they still PASS with a wrong V_G?
(B) NM-3: independent a(s) by quadrature of G vs the closed-form edge law, and domain sensitivity.
(C) reviewer's own physical-time residual for the pure-step equation with real-space Q and V_G by quadrature."""
import sys, io, math, contextlib, numpy as np
np.seterr(all="ignore")
REPO = "/Volumes/4TB SD/ClaudeCode/singular-blowup-1d/validate/"
import wk, rsq
which = sys.argv[1]
def run_src(src, fn):
    g = {"__name__": "gate_copy"}
    exec(compile(src, "gate_copy", "exec"), g)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        try: ok = g[fn]()
        except Exception as e: ok = f"EXC {e!r}"
    return ok, buf.getvalue().strip(), g
if which == "A":
    pm = open(REPO + "edge_pure_step_nonlinear_manifold.py").read()
    nm = open(REPO + "edge_nonlinear_cusp_manifold.py").read()
    pm_ret = "return float(-(1.0 / (math.pi * r_val)) * np.sum(np.exp(0.5 * z) * ln_term * G_prof) * dz)"
    nm_ret = "return float(-(1.0 / (math.pi * r_val)) * np.sum(ig) * dz)"
    assert pm_ret in pm and nm_ret in nm
    muts = {"unmutated": "{R}", "V_G -> 3*V_G": "{R3}", "V_G -> -V_G": "{RM}", "V_G -> 0.05*sin(7*r) (unrelated to G)": "return 0.05*math.sin(7.0*r_val)",
            "V_G -> 0.3*r": "return 0.3*r_val"}
    for name, m in muts.items():
        for label, src, ret, fn in (("PM-2/3", pm, pm_ret, "check_pm2_pm3"), ("NM-3", nm, nm_ret, "check_nm3")):
            new = m.replace("{R3}", ret.replace("return float(", "return 3.0*float(")).replace("{RM}", ret.replace("return float(", "return -1.0*float(")).replace("{R}", ret)
            ok, out, _ = run_src(src.replace(ret, new), fn)
            print(f"[{name}] {label}: returns {ok}\n   " + out.replace("\n", "\n   "))
if which == "B":
    nm = open(REPO + "edge_nonlinear_cusp_manifold.py").read()
    hook = "    ok = bool(\n        diff_G < 1e-10"
    assert hook in nm
    for (zlo, zhi, Nz) in ((-32, 32, 2048), (-48, 48, 3072), (-64, 64, 4096), (-96, 32, 4096)):
        src = nm.replace(hook, '    globals()["_L"] = dict(locals())\n' + hook).replace("ops = build_spectral_ops(Nz=Nz, zmin=-32.0, zmax=32.0)", f"ops = build_spectral_ops(Nz=Nz, zmin={zlo}.0, zmax={zhi}.0)").replace("    Nz = 2048\n    ops", f"    Nz = {Nz}\n    ops")
        ok, out, g = run_src(src, "check_nm3")
        print(f"--- NM-3 on z=[{zlo},{zhi}] Nz={Nz}: returns {ok}\n   " + out.replace("\n", "\n   "))
        Lc = g.get("_L")
        if Lc is None: continue
        z, G, H, QG, s, W = Lc["z"], Lc["G"], Lc["H"], Lc["QG_all"], Lc["s_grid"], Lc["W"]
        a0 = -1 / (2 * 2**1.5) * 0.05
        print("   s : a_closed=a0/(1+a0 s/2) | a_quad(G)[-30,20] | Qconj(G) at z=-10,-14,-20,-26 | -H at z=-14,-20 | W at z=-14,-20 | max|G| z in(12,zhi-2)")
        for j in (0, 1, 8, 48, 96):
            ac = a0 / (1 + 0.5 * a0 * s[j])
            print(f"   {s[j]:.4f}: {ac:+.8e} | {rsq.a_quad(G[j], z, -30, 20):+.8e} | " + " ".join(f"{QG[j][np.argmin(np.abs(z-zz))]:+.8e}" for zz in (-10, -14, -20, -26))
                  + " | " + " ".join(f"{-H[j][np.argmin(np.abs(z-zz))]:+.8e}" for zz in (-14, -20)) + " | " + " ".join(f"{W[j][np.argmin(np.abs(z-zz))]:+.2e}" for zz in (-14, -20))
                  + f" | {np.max(np.abs(G[j][(z>12)&(z<zhi-2)])):.2e}")
if which == "C":
    ops = wk.Ops(Nz=8192, zmin=-80, zmax=80, filt_frac=None); z, dz = ops.z, ops.dz
    Ns = 480; ds = 1.2 / Ns; H = {}
    def rec(j, s, h): H[j] = ops.Ng(h)
    wk.etdrk2(ops, wk.seed_u3(ops), eps=0.02, Ns=Ns, rec_every=1, rec=rec)
    m = (z >= -40) & (z <= 30); zz = z[m]
    def V_of(Gt, r):      # V_G of the PHYSICAL profile G(v)=Gt(r^2 v):  -(1/pi) int v^-1/2 ln(1+1/v) G dv, in w=r^2 v variables
        return -(1 / (math.pi * r)) * np.sum(np.exp(0.5 * zz) * (np.log1p(r * r * np.exp(-zz))) * Gt[m]) * dz
    r = np.zeros(Ns + 1); r[0] = 1.0; V = np.zeros(Ns + 1); V[0] = V_of(H[0], 1.0); jend = Ns
    for j in range(Ns):             # Heun
        k1 = -1 + 0.5 * V_of(H[j], r[j]); rp = r[j] + ds * k1
        if rp <= 1e-3: jend = j; break
        k2 = -1 + 0.5 * V_of(H[j + 1], rp); r[j + 1] = r[j] + 0.5 * ds * (k1 + k2); V[j + 1] = V_of(H[j + 1], r[j + 1])
        if r[j + 1] <= 1e-3: jend = j; break
    tau = np.concatenate([[0], np.cumsum(0.5 * ds * (1 / r[:jend] + 1 / r[1:jend + 1]))])
    print(f"C: r hits ~0 near s={jend*ds:.4f} (gate s_inf=0.998846); tau reached {tau[-1]:.2f}")
    from scipy.interpolate import CubicSpline
    zv = np.arange(-14.0, 10.0, dz * 2); Hm = rsq.H_matrix(np.arange(-34.0, 26.0, dz * 2)); zq = np.arange(-34.0, 26.0, dz * 2)
    for tt in (0.2, 0.6, 1.5, 3.0):
        j = int(np.argmin(np.abs(tau - tt)))
        def phys(jj, grid): return CubicSpline(z, H[jj])(grid + 2 * math.log(r[jj]))
        Gq = phys(j, zq); q = rsq.EA_apply(Hm @ Gq, zq)            # real-space Q of the physical profile
        Vq = -(1 / math.pi) * np.sum(np.exp(0.5 * zq) * np.log1p(np.exp(-zq)) * Gq) * (2 * dz)   # V_G directly in physical v
        dtau = tau[j + 1] - tau[j - 1]
        Gt = (phys(j + 1, zq) - phys(j - 1, zq)) / dtau
        # non-uniform tau spacing correction: second-order centred at tau_j only if symmetric; use 3-pt Lagrange derivative
        t0, t1, t2 = tau[j - 1], tau[j], tau[j + 1]
        Gt = (phys(j - 1, zq) * (t1 - t2) / ((t0 - t1) * (t0 - t2)) + Gq * (2 * t1 - t0 - t2) / ((t1 - t0) * (t1 - t2)) + phys(j + 1, zq) * (t1 - t0) / ((t2 - t0) * (t2 - t1)))
        Gz = rsq.D6(Gq, 2 * dz)
        res = Gt + (2 - Vq + q) * Gz - 0.5 * q * (1 + Gq)
        lin = Gt + 2 * Gz - 0.5 * q
        for nm_, mk in (("edge v in[e^-14,e^-8]", (zq >= -14) & (zq < -8)), ("bulk [e^-8,e^8]", (zq >= -8) & (zq <= 8)), ("far (e^8,e^12]", (zq > 8) & (zq <= 12))):
            print(f"   tau={tau[j]:.3f} (s={j*ds:.4f}, r={r[j]:.4f}) V_G={Vq:+.5e} (w-form {V[j]:+.5e}) {nm_:22s}: max|res|={np.nanmax(np.abs(res[mk])):.2e}  [dropping nonlinear+V terms: {np.nanmax(np.abs(lin[mk])):.2e}]  max|G|={np.max(np.abs(Gq[mk])):.2e}  Q(edge)={q[0]:+.1e}")
