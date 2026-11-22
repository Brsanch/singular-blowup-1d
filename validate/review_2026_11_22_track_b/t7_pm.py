from mut import run
VG_OLD="        return float(-(1.0 / (math.pi * r_val)) * np.sum(np.exp(0.5 * z) * ln_term * G_prof) * dz)\n\n    r_arr"
def vg(new): return [(VG_OLD, "        return "+new+"\n\n    r_arr")]
BASE="float(-(1.0 / (math.pi * r_val)) * np.sum(np.exp(0.5 * z) * ln_term * G_prof) * dz)"
dom=lambda zmin,zmax,nz:[("ops = build_ops(Nz=4096, zmin=-64.0, zmax=64.0)\n    z, w, u, dz = ops",f"ops = build_ops(Nz={nz}, zmin={zmin}, zmax={zmax})\n    z, w, u, dz = ops"),("np.fft.fftfreq(4096, dz)",f"np.fft.fftfreq({nz}, dz)")]
run("PM-2/3 baseline","PM","check_pm2_pm3",[])
print("\n##### C3: consistent V_G mutations (formula changed everywhere it is used)")
run("PM-3 MUT V_G -> 3 V_G (consistently)","PM","check_pm2_pm3",vg("3.0*"+BASE))
run("PM-3 MUT V_G -> -V_G (consistently)","PM","check_pm2_pm3",vg("-1.0*"+BASE))
run("PM-3 MUT V_G -> 0.05 sin(7 r) (unrelated function)","PM","check_pm2_pm3",vg("0.05*math.sin(7.0*r_val)"))
run("PM-3 MUT V_G -> 0 (consistently)","PM","check_pm2_pm3",vg("0.0*"+BASE))
print("\n##### C1: a-metric under wrong nonlinearity in the solver (and residual, consistently)")
NLs=("fP[j] = NgI(-QGj * Dg(Gj) + 0.5 * QGj * Gj)","float(np.max(np.abs(((G[j + 1] - G[j - 1]) / (2.0 * ds) - 0.5 * QG_all[j] + QG_all[j] * Dg(G[j]) - 0.5 * QG_all[j] * G[j])[bulk])))")
run("PM-2 MUT advection sign flipped (solver + residual)","PM","check_pm2_pm3",[(NLs[0],"fP[j] = NgI(+QGj * Dg(Gj) + 0.5 * QGj * Gj)"),(NLs[1],NLs[1].replace("+ QG_all[j] * Dg(G[j])","- QG_all[j] * Dg(G[j])"))])
run("PM-2 MUT nonlinearity -> G^2 (no Q) (solver + residual)","PM","check_pm2_pm3",[(NLs[0],"fP[j] = NgI(Gj * Gj)"),(NLs[1],NLs[1].replace("+ QG_all[j] * Dg(G[j]) - 0.5 * QG_all[j] * G[j]","- G[j] * G[j]"))])
run("PM-2 MUT advection sign flipped in solver ONLY (residual = true equation)","PM","check_pm2_pm3",[(NLs[0],"fP[j] = NgI(+QGj * Dg(Gj) + 0.5 * QGj * Gj)")])
run("PM-2 MUT wrong Q operator Q = -2 Nr(u^1.5 NgI) everywhere","PM","check_pm2_pm3",[("    Q = lambda G: -2.0 * Nr(u * NgI(G))\n    Dg = lambda G: np.fft.ifft(m_g * np.fft.fft(w_g * G)).real / w_g","    Q = lambda G: -2.0 * Nr(u**1.5 * NgI(G))\n    Dg = lambda G: np.fft.ifft(m_g * np.fft.fft(w_g * G)).real / w_g")])
print("\n##### C1/C2: N_s, N_z, domain convergence of the printed numbers")
run("PM Ns=480","PM","check_pm2_pm3",[("    Ns = 240\n    s_grid = np.linspace(0.0, S_max, Ns + 1)","    Ns = 480\n    s_grid = np.linspace(0.0, S_max, Ns + 1)"),("abs(a_quad_tail(G[120])), abs(a_quad_tail(G[200]))","abs(a_quad_tail(G[240])), abs(a_quad_tail(G[400]))"),("for j in range(10, Ns)","for j in range(20, Ns)")])
run("PM Ns=960","PM","check_pm2_pm3",[("    Ns = 240\n    s_grid = np.linspace(0.0, S_max, Ns + 1)","    Ns = 960\n    s_grid = np.linspace(0.0, S_max, Ns + 1)"),("abs(a_quad_tail(G[120])), abs(a_quad_tail(G[200]))","abs(a_quad_tail(G[480])), abs(a_quad_tail(G[800]))"),("for j in range(10, Ns)","for j in range(40, Ns)")])
run("PM Nz=8192 on [-64,64]","PM","check_pm2_pm3",dom(-64.0,64.0,8192))
run("PM Nz=2048 on [-64,64]","PM","check_pm2_pm3",dom(-64.0,64.0,2048))
run("PM domain [-48,48] Nz=4096","PM","check_pm2_pm3",dom(-48.0,48.0,4096))
run("PM domain [-96,96] Nz=8192","PM","check_pm2_pm3",dom(-96.0,96.0,8192))
run("PM domain [-96,96] Nz=4096","PM","check_pm2_pm3",dom(-96.0,96.0,4096))
run("PM domain [-104,24] Nz=4096 (long-left, C5)","PM","check_pm2_pm3",dom(-104.0,24.0,4096))
run("PM eps=0.02 (the value the docstring/note quote)","PM","check_pm2_pm3",[("    eps = 0.05\n    S_max = 1.2","    eps = 0.02\n    S_max = 1.2")])
run("PM eps=0.2","PM","check_pm2_pm3",[("    eps = 0.05\n    S_max = 1.2","    eps = 0.2\n    S_max = 1.2")])
