import sys, numpy as np
from scipy.linalg import expm
sys.path.insert(0, '/Volumes/4TB SD/ClaudeCode/singular-blowup-1d/src')
from stability_edge_z import build_blocks_z
GAMMA, Z0, S = 0.4, -6.0, 2.0
TAUS = np.linspace(0.0, 6.0, 31)
print("transient_growth_edge_conditioned.py, reduced reproduction (same formulas), tail off vs on")
for (N, zmin, zmax) in ((800, -10, 8), (1000, -12, 8)):
    for tail in (False, True):
        B = build_blocks_z(N=N, zmin=zmin, zmax=zmax, beta=0.0, edge_bc=True, edge_tail=tail)
        M = -2.0*B["Dz"] + 0.5*B["Q_op"]
        wc = max(float(np.linalg.norm(expm(M*t), 2)) for t in TAUS)
        z, dz, Qop, Dz = B["z"], B["dz"], B["Q_op"], B["Dz"]
        r = []
        for z_cut in (-8.0, -6.0):
            keep = z > z_cut; zk = z[keep]
            pp = 2.0*GAMMA*0.5*(1.0+np.tanh((zk-Z0)/S)); sw = np.exp(-0.5*(np.cumsum(pp)*dz))
            Fb = 0.5*Qop[np.ix_(keep, keep)]; Ft = (sw[:, None]*Fb)/sw[None, :]
            mu_an = float(np.linalg.eigvalsh(-np.diag(pp)+0.5*(Ft+Ft.T))[-1])
            Mb = -2.0*Dz[np.ix_(keep, keep)] + Fb; Mt = (sw[:, None]*Mb)/sw[None, :]
            mu_disc = float(np.linalg.eigvalsh(0.5*(Mt+Mt.T))[-1])
            r.append(f"z_cut={z_cut:+.0f}: mu_an={mu_an:+.5f} mu_disc={mu_disc:+.2f}")
        print(f" N={N} zmin={zmin} tail={tail!s:5}: worst-case sup||e^(Mt)||_2={wc:8.2f} | " + " ; ".join(r)); sys.stdout.flush()
