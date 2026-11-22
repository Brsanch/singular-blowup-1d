#!/usr/bin/env python3
"""Stability / domain / filter scan of the gates' conjugated ETD scheme (PM seed u^3 e^-u and FS bump).
Reports blow-up time, far-field garbage, edge slope, a by quadrature (tail-corrected)."""
import sys, math, numpy as np, wk, rsq
np.seterr(all="ignore")
dz0 = 80.0 / 4096
def run(seed, eps, ff, zlo, zhi, Ns=240, nl_name="P", dzf=1.0):
    Nz = int(round((zhi - zlo) / (dz0 * dzf)))
    Nz += Nz % 2
    ops = wk.Ops(Nz=Nz, zmin=zlo, zmax=zhi, filt_frac=ff)
    z, u, w = ops.z, ops.u, ops.w
    h0 = wk.seed_u3(ops) if seed == "u3" else wk.seed_bump(ops)
    nl = {"P": ops.P, "G2": lambda G: 5.0 * G * G, "bumpG": lambda G: np.exp(-0.5 * (z + 3) ** 2) * G * 0.5,
          "linP": None}[nl_name]
    out = {}
    blow = [None]
    def rec(j, s, h):
        if blow[0] is None and (not np.all(np.isfinite(h)) or np.max(np.abs(h[(z > -12) & (z < 12)])) > 1e3):
            blow[0] = s
        out[round(s, 3)] = h.copy()
    h = wk.etdrk2(ops, h0, eps=eps, Ns=Ns, nl=nl, rec_every=max(1, Ns // 24), rec=rec)
    tag = f"{seed} eps={eps} nl={nl_name} filt={ff} z=[{zlo},{zhi}] Nz={Nz} Ns={Ns}"
    if blow[0] is not None:
        print(f"{tag}: BLOW-UP/NaN first seen at s={blow[0]:.2f}")
        return
    for sv in (0.6, 1.2):
        hh = out[sv]; G = ops.Ng(hh)
        core = (z > -8) & (z < 8)
        far = (z > 12) & (z < zhi - 2)
        aP = rsq.a_quad(nl(G), z, -20, min(20, zhi - 2))
        k14 = np.argmin(np.abs(z + 14)); k20 = np.argmin(np.abs(z + 20))
        c1 = G[k14] / w[k14]
        aq = {lo: rsq.a_quad(G, z, lo, min(20, zhi - 2)) for lo in (-20, -30)}
        aq_corr = {lo: aq[lo] - (1 / math.pi) * c1 * 2 * math.exp(lo / 2) for lo in aq}   # add analytic c1*w tail below lo
        ie = np.argmin(np.abs(z + 15))
        print(f"{tag} s={sv}: max|G|core={np.max(np.abs(G[core])):.6e} max|G|far(z>12)={np.max(np.abs(G[far])):.2e} "
              f"G/w(-14)={c1:+.5e} G/w(-20)={G[k20]/w[k20]:+.5e} 3a(nl)={3*aP:+.5e} | gate a=Qconj(-15)={ops.Q(G)[ie]:+.2e} "
              f"a_quad[-20]={aq[-20]:+.2e} a_quad[-30]={aq[-30]:+.2e} tail-corrected: {aq_corr[-20]:+.2e} {aq_corr[-30]:+.2e} "
              f"| norms L2,H1,H2,H3 = " + " ".join(f"{x:.5e}" for x in ops.norms(hh)))
if __name__ == "__main__":
    which = sys.argv[1]
    if which == "pm":
        for ff in (0.6, None):
            for (lo, hi) in ((-40, 40), (-48, 48), (-56, 56), (-64, 64), (-80, 80), (-104, 24), (-64, 24)):
                run("u3", 0.02, ff, lo, hi)
        run("u3", 0.002, 0.6, -104, 24)
        run("u3", 0.02, 0.6, -80, 80, nl_name="G2")
    if which == "mut":
        for nlm in ("P", "G2", "bumpG"):
            run("u3", 0.02, 0.6, -40, 40, nl_name=nlm)
            run("u3", 0.02, None, -64, 24, nl_name=nlm)
