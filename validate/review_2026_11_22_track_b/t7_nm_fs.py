from mut import run, SRC
print("##### PM: does an end-sponge on the forcing cure the Nz=8192 NaN?  (mechanism test on the gate's own Picard solver)")
dom=lambda zmin,zmax,nz:[("ops = build_ops(Nz=4096, zmin=-64.0, zmax=64.0)\n    z, w, u, dz = ops",f"ops = build_ops(Nz={nz}, zmin={zmin}, zmax={zmax})\n    z, w, u, dz = ops"),("np.fft.fftfreq(4096, dz)",f"np.fft.fftfreq({nz}, dz)")]
SP=[("            fP[j] = NgI(-QGj * Dg(Gj) + 0.5 * QGj * Gj)\n","            fP[j] = NgI(-QGj * Dg(Gj) + 0.5 * QGj * Gj) * (0.5*(1+np.tanh((z-z[0]-8.0)/1.5))) * (0.5*(1+np.tanh((z[-1]-8.0-z)/1.5)))\n")]
run("PM Nz=8192 [-64,64] + tanh end-sponge on fP","PM","check_pm2_pm3",dom(-64.0,64.0,8192)+SP)
run("PM [-104,24] Nz=4096 + tanh end-sponge on fP","PM","check_pm2_pm3",dom(-104.0,24.0,4096)+SP)

print("\n##### NM-3")
VG_OLD="        ig = np.exp(0.5 * z) * ln_term * G_prof\n        return float(-(1.0 / (math.pi * r_val)) * np.sum(ig) * dz)"
def vg(new): return [(VG_OLD, "        ig = np.exp(0.5 * z) * ln_term * G_prof\n        return "+new)]
B="float(-(1.0 / (math.pi * r_val)) * np.sum(ig) * dz)"
run("NM-3 baseline","NM","check_nm3",[])
run("NM-3 MUT V_G -> 3 V_G (consistently)","NM","check_nm3",vg("3.0*"+B))
run("NM-3 MUT V_G -> -V_G (consistently)","NM","check_nm3",vg("-1.0*"+B))
run("NM-3 MUT V_G -> 0.05 sin(7 r)","NM","check_nm3",vg("0.05*math.sin(7.0*r_val)"))
nmdom=lambda zmin,zmax,nz:[("    Nz = 4096\n    ops = build_spectral_ops(Nz=Nz, zmin=-64.0, zmax=64.0)",f"    Nz = {nz}\n    ops = build_spectral_ops(Nz=Nz, zmin={zmin}, zmax={zmax})")]
run("NM-3 Nz=8192 [-64,64]","NM","check_nm3",nmdom(-64.0,64.0,8192))
run("NM-3 [-96,96] Nz=8192","NM","check_nm3",nmdom(-96.0,96.0,8192))
run("NM-3 [-48,48] Nz=4096","NM","check_nm3",nmdom(-48.0,48.0,4096))
run("NM-3 [-104,24] Nz=4096","NM","check_nm3",nmdom(-104.0,24.0,4096))
run("NM-3 Ns=256","NM","check_nm3",[("    Ns = 128\n    Nz = 4096","    Ns = 256\n    Nz = 4096")])
run("NM-3 filter off (0.6 -> 1.0)","NM","check_nm3",[("filt = np.where(np.abs(xi) <= 0.6 * xi_max, 1.0, 0.0)","filt = np.where(np.abs(xi) <= 1.0 * xi_max, 1.0, 0.0)")])
run("NM-3 eps=0.05 (value quoted in docstring/note)","NM","check_nm3",[("    eps = 0.08\n    S_max = 1.2\n    Ns = 128","    eps = 0.05\n    S_max = 1.2\n    Ns = 128")])
run("NM-2 baseline (still [-36,36], filter on)","NM","check_nm2",[])
run("NM-2 on [-64,64]","NM","check_nm2",[("ops = build_spectral_ops(Nz=4096, zmin=-36.0, zmax=36.0)","ops = build_spectral_ops(Nz=4096, zmin=-64.0, zmax=64.0)")])
run("NM-2 MUT wrong source N2 -> 3*N2: does the 'automatic edge balance' still pass?","NM","check_nm2",[("    N2 = 0.5 * ops[\"half_plus_D_r\"](q0_exact**2) + 0.5 * QS1\n","    N2 = 3.0*(0.5 * ops[\"half_plus_D_r\"](q0_exact**2) + 0.5 * QS1) + 0.1*np.exp(-(z+2)**2)\n")])

print("\n##### FS-3")
fsdom=lambda zmin,zmax:[("        zmin, zmax = -72.0, 72.0\n",f"        zmin, zmax = {zmin}, {zmax}\n")]
run("FS-3 baseline","FS","check_fs3",[])
run("FS-3 Nz in (8192,16384) on [-72,72]","FS","check_fs3",[("for Nz in (4096, 8192):","for Nz in (8192, 16384):"),("results[8192][i] - results[4096][i]) / results[8192][i]","results[16384][i] - results[8192][i]) / results[16384][i]"),("a_q = abs(results[8192][4])","a_q = abs(results[16384][4])"),("results[8192]","results[16384]")])
run("FS-3 gate-original domain [-36,36]","FS","check_fs3",fsdom(-36.0,36.0))
run("FS-3 [-64,64]","FS","check_fs3",fsdom(-64.0,64.0))
run("FS-3 [-96,96]","FS","check_fs3",fsdom(-96.0,96.0))
run("FS-3 [-104,24] (long-left)","FS","check_fs3",fsdom(-104.0,24.0))
run("FS-3 Ns=480","FS","check_fs3",[("    Ns = 240\n    ds = S_max / Ns\n\n    results","    Ns = 480\n    ds = S_max / Ns\n\n    results")])
run("FS-3 eps=0.06","FS","check_fs3",[("    eps = 0.03\n    S_max = 1.2\n    Ns = 240\n    ds = S_max / Ns\n\n    results","    eps = 0.06\n    S_max = 1.2\n    Ns = 240\n    ds = S_max / Ns\n\n    results")])
run("FS-3 eps=0.12","FS","check_fs3",[("    eps = 0.03\n    S_max = 1.2\n    Ns = 240\n    ds = S_max / Ns\n\n    results","    eps = 0.12\n    S_max = 1.2\n    Ns = 240\n    ds = S_max / Ns\n\n    results")])
run("FS-3 MUT nonlinearity -> G^2: a_quad still ~0?","FS","check_fs3",[("            return NgI(-QG_cur * Dg(G_cur) + 0.5 * QG_cur * G_cur)","            return NgI(G_cur * G_cur)")])
run("FS-1 baseline","FS","check_fs1",[])
run("FS-1 with k=0,1 restored (A4)","FS","check_fs1",[("for m_ord in (2, 3):","for m_ord in (0, 1, 2, 3):")])
