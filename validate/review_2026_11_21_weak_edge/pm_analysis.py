#!/usr/bin/env python3
"""Tasks 3+4 on the gate's own PM-2 solution: independent a(s) and real-space residual."""
import sys, math
import numpy as np
import wk, rsq

np.set_printoptions(linewidth=200)
ff = float(sys.argv[1]) if len(sys.argv) > 1 else 0.6
zlo_ = float(sys.argv[2]) if len(sys.argv) > 2 else -40.0
zm = float(sys.argv[3]) if len(sys.argv) > 3 else 40.0
Nz = int(sys.argv[4]) if len(sys.argv) > 4 else 4096
seed = sys.argv[5] if len(sys.argv) > 5 else "u3"
eps = float(sys.argv[6]) if len(sys.argv) > 6 else 0.02
ops = wk.Ops(Nz=Nz, zmin=zlo_, zmax=zm, filt_frac=(None if ff <= 0 else ff))
z, u, dz = ops.z, ops.u, ops.dz
h0 = wk.seed_u3(ops) if seed == "u3" else wk.seed_bump(ops)
Ns = 240
s, G, hh, it, diff = wk.picard(ops, h0, eps=eps, Ns=Ns)
ds = s[1] - s[0]
bulk = (z > -12) & (z < 12)
QG = np.array([ops.Q(g) for g in G])
res_gate = max(float(np.max(np.abs(((G[j+1]-G[j-1])/(2*ds) - 0.5*QG[j] + QG[j]*ops.Dg(G[j]) - 0.5*QG[j]*G[j])[bulk]))) for j in range(10, Ns))
ie = np.argmin(np.abs(z + 15.0))
print(f"# run filt={ff} zmin={zlo_} zmax={zm} Nz={Nz} seed={seed} eps={eps}: iter={it} diff={diff:.2e} gate-res={res_gate:.3e} gate max|a|={np.max(np.abs(QG[:, ie])):.3e}")
jm = int(np.argmax(np.abs(QG[:, ie])))
print(f"#   gate max|a| attained at s={s[jm]:.3f}; excluding s=0: {np.max(np.abs(QG[1:, ie])):.3e}")

# --- sizes of terms (how discriminating is a 1e-4 threshold?)
j = 120
Pj = ops.P(G[j])
print(f"# term sizes at s={s[j]:.2f} (bulk max): |G|={np.max(np.abs(G[j][bulk])):.3e} |0.5QG|={0.5*np.max(np.abs(QG[j][bulk])):.3e} "
      f"|QG*DG|={np.max(np.abs((QG[j]*ops.Dg(G[j]))[bulk])):.3e} |0.5QG*G|={0.5*np.max(np.abs((QG[j]*G[j])[bulk])):.3e} |P|={np.max(np.abs(Pj[bulk])):.3e}")

# --- mutation of the residual evaluator (gate metric, same solution, WRONG pde)
def gate_res(c_adv, c_quad, c_lin=0.5, j0=10):
    return max(float(np.max(np.abs(((G[j+1]-G[j-1])/(2*ds) - c_lin*QG[j] + c_adv*QG[j]*ops.Dg(G[j]) - c_quad*QG[j]*G[j])[bulk]))) for j in range(j0, Ns))
print("# gate-metric residual of the SAME solution against mutated PDEs (threshold 1e-4):")
for ca, cq in ((1, 0.5), (1, 0.4), (1, 0.0), (0.9, 0.5), (0.8, 0.5), (0.0, 0.5), (0.0, 0.0)):
    print(f"#   c_adv={ca} c_quad={cq}: {gate_res(ca, cq):.3e}")
# where/when is the gate residual attained
best = (0, 0, 0)
for j in range(10, Ns):
    r = np.abs((G[j+1]-G[j-1])/(2*ds) - 0.5*QG[j] + QG[j]*ops.Dg(G[j]) - 0.5*QG[j]*G[j]) * bulk
    k = int(np.argmax(r))
    if r[k] > best[0]:
        best = (r[k], s[j], z[k])
print(f"#   gate residual max {best[0]:.3e} at s={best[1]:.3f}, z={best[2]:.2f}")
# exact-in-s version (no FD in s): G_s = Ng(-u h + f)
def res_exact_s(j, Qv=None, Dv=None, mask=bulk):
    f = ops.NgI(ops.P(G[j]))
    Gs = ops.Ng(-u * hh[j] + f)
    Qv = QG[j] if Qv is None else Qv
    Dv = ops.Dg(G[j]) if Dv is None else Dv
    return Gs - 0.5*Qv + Qv*Dv - 0.5*Qv*G[j]
print(f"#   same residual with exact G_s=Ng(-u h+f) (pure algebraic identity, expected ~roundoff): "
      f"{max(float(np.max(np.abs(res_exact_s(j)[bulk]))) for j in range(10, Ns)):.3e}")

# --- Task 3: a three ways
print("\n# TASK 3: a(s) three ways.  conj = Q_conj(G)(z0);  quad = -(1/pi) sum e^{-z/2} G dz on [lo,hi];  rsQ = real-space Q edge value")
wins = [(-20, 20), (-25, 25), (-30, 30), (-36, 36)]
wins = [(lo, min(hi, zm - 2)) for lo, hi in wins]
Hms = {}
for (lo, hi) in [(-20, min(20, zm - 2)), (-28, min(28, zm - 2))]:
    m = (z >= lo) & (z <= hi)
    Hms[(lo, hi)] = (m, rsq.H_matrix(z[m]))
hdr = "s      " + " ".join(f"conj(z={zz:+d})" for zz in (-10, -15, -20, -25, -30)) + " | " + " ".join(f"quad[{lo},{hi}]" for lo, hi in wins) + " | quad[full] | " + " ".join(f"rsQ[{lo},{hi}](edge) rsQ(z=-15)" for lo, hi in Hms)
print(hdr)
for j in (0, 1, 2, 10, 40, 120, 200, 240):
    row = f"{s[j]:.3f} "
    for zz in (-10, -15, -20, -25, -30):
        row += f" {QG[j, np.argmin(np.abs(z - zz))]:+.3e}"
    row += " |"
    for lo, hi in wins:
        row += f" {rsq.a_quad(G[j], z, lo, hi):+.3e}"
    row += f" | {rsq.a_quad(G[j], z):+.3e} |"
    for key, (m, Hm) in Hms.items():
        q = rsq.EA_apply(Hm @ G[j][m], z[m])
        row += f" {q[0]:+.3e} {q[np.argmin(np.abs(z[m] + 15.0))]:+.3e}"
    print(row)
# scale references
print(f"# reference scales: (1/pi) int e^(-z/2)|G| dz at s=0.6: {np.sum(np.exp(-0.5*z[(z>-30)&(z<30)])*np.abs(G[120][(z>-30)&(z<30)]))*dz/np.pi:.3e};  a(P(G)) at s=0.6 by quad[-20,20]: {rsq.a_quad(ops.P(G[120]), z, -20, 20):+.3e}")
# edge slope c1 vs 3 a(P)
for j in (40, 120, 240):
    for zz in (-8, -10, -12, -14):
        k = np.argmin(np.abs(z - zz))
        print(f"#   s={s[j]:.2f} z={zz}: G/w={G[j][k]/ops.w[k]:+.6e}   3*a(P) quad[-20,20]={3*rsq.a_quad(ops.P(G[j]), z, -20, 20):+.6e}")

# --- Task 4: real-space residual
print("\n# TASK 4: real-space residual  G_s - Q/2 + Q*DG - Q*G/2,  Q = reviewer p.v. quadrature, D = 6th-order FD")
for key, (m, Hm) in Hms.items():
    zz = z[m]
    reg = {"edge[-18,-12]": (zz >= -18) & (zz < -12), "bulk[-12,12]": (zz >= -12) & (zz <= 12), "far(12,18]": (zz > 12) & (zz <= 18), "core[-8,8]": (zz >= -8) & (zz <= 8)}
    out = {k: [0.0, 0.0, 0.0] for k in reg}
    dq = {k: 0.0 for k in reg}
    for j in list(range(10, Ns, 10)):
        Gj = G[j][m]
        q = rsq.EA_apply(Hm @ Gj, zz)
        D = rsq.D6(Gj, dz)
        Gs_fd = (G[j+1][m] - G[j-1][m]) / (2*ds)
        Gs_ex = ops.Ng(-u * hh[j] + ops.NgI(ops.P(G[j])))[m]
        r_fd = Gs_fd - 0.5*q + q*D - 0.5*q*Gj
        r_ex = Gs_ex - 0.5*q + q*D - 0.5*q*Gj
        r_lin = Gs_fd - 0.5*q            # nonlinear terms dropped: what a WRONG pde gives
        for k, mk in reg.items():
            out[k][0] = max(out[k][0], float(np.nanmax(np.abs(r_fd[mk]))))
            out[k][1] = max(out[k][1], float(np.nanmax(np.abs(r_ex[mk]))))
            out[k][2] = max(out[k][2], float(np.nanmax(np.abs(r_lin[mk]))))
            dq[k] = max(dq[k], float(np.max(np.abs((q - QG[j][m])[mk]))))
    for k in reg:
        print(f"window {key} {k:14s}: res(FD-in-s)={out[k][0]:.3e}  res(exact G_s)={out[k][1]:.3e}  [linear-only pde: {out[k][2]:.3e}]  max|Q_rs-Q_conj|={dq[k]:.3e}")
np.savez("/private/tmp/claude-501/-Volumes-4TB-SD-ClaudeCode-noethersolve/24da0ae4-02c7-4541-b503-95490895dcb6/scratchpad/weakedge_review/pm_sol_%s_%g_%g_%g_%d.npz" % (seed, ff, zlo_, zm, Nz), s=s, z=z, G=G, h=hh)

print("# profile at s=0.6: z, G, G/w, h*u^2")
for zz in (-30,-24,-20,-16,-14,-12,-10,-8,10,12,14,16,18,22):
    if zz < zm-1:
        k=np.argmin(np.abs(z-zz)); print(f"#   z={zz:+d} G={G[120][k]:+.4e} G/w={G[120][k]/ops.w[k]:+.5e} u^2h={u[k]**2*hh[120][k]:+.5e}")
print(f"#   targets at s=0.6: 3a(P)={3*rsq.a_quad(ops.P(G[120]), z, -20, min(20,zm-2)):+.5e}  -pi/2 a(P)={-0.5*math.pi*rsq.a_quad(ops.P(G[120]), z, -20, min(20,zm-2)):+.5e}")
