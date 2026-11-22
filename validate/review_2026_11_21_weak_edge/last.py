import sys, io, contextlib, numpy as np, wk, rsq
np.seterr(all="ignore")
REPO = "/Volumes/4TB SD/ClaudeCode/singular-blowup-1d/validate/"
nm = open(REPO + "edge_nonlinear_cusp_manifold.py").read()
hook = "    ok = bool(\n        diff_G < 1e-10"
for Ns in (48, 96, 192, 384):
    src = nm.replace(hook, '    globals()["_L"] = dict(locals())\n' + hook).replace("ops = build_spectral_ops(Nz=Nz, zmin=-32.0, zmax=32.0)", "ops = build_spectral_ops(Nz=Nz, zmin=-64.0, zmax=64.0)").replace("    Nz = 2048\n    ops", "    Nz = 4096\n    ops").replace("    Ns = 96\n", f"    Ns = {Ns}\n")
    g = {"__name__": "gc"}; exec(compile(src, "gc", "exec"), g)
    with contextlib.redirect_stdout(io.StringIO()) as b: ok = g["check_nm3"]()
    L = g["_L"]; z = L["z"]; j = L["Ns"] // 2
    aq = rsq.a_quad(L["G"][j], z, -30, 20); D = L["H"][j][np.argmin(np.abs(z + 20))]
    print(f"NM-3 clean domain Ns={Ns}: gate ok={ok} |a+D| gate={L['max_s_compat']:.3e}; s=0.6: a_quad(G)+D(z=-20) = {aq + D:+.3e}; res_G={L['max_s_res_G']:.2e}")
# PM: early-time gate 'a' vs time step (clean domain, Picard-trapezoid as in gate)
for Ns in (240, 480, 960, 1920):
    ops = wk.Ops(Nz=8192, zmin=-80, zmax=80, filt_frac=None); z = ops.z
    s, G, hh, it, diff = wk.picard(ops, wk.seed_u3(ops), eps=0.02, Ns=Ns)
    QG = np.array([ops.Q(g) for g in G[: Ns // 4 + 1]])
    row = []
    for zz in (-12, -15, -18):
        k = np.argmin(np.abs(z - zz)); jm = int(np.argmax(np.abs(QG[1:, k]))) + 1
        row.append(f"z={zz}: max|Qconj|={abs(QG[jm, k]):.2e} at s={s[jm]:.4f} (u*ds={ops.u[k]*1.2/Ns:.1f})")
    print(f"PM clean domain Ns={Ns} iters={it}: " + "; ".join(row))
