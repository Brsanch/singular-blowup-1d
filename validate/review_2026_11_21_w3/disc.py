# conjugation-free: full weighted operator norm of the repo's discrete semigroup
import sys, numpy as np
from scipy.linalg import expm
sys.path.insert(0, "/Volumes/4TB SD/ClaudeCode/singular-blowup-1d/src")
from stability_edge_z import build_blocks_z
for (N,zmin,zmax) in [(1000,-10.,10.),(1400,-14.,10.),(1400,-18.,10.)]:
    B=build_blocks_z(N=N,zmin=zmin,zmax=zmax,edge_bc=True); z,M=B["z"],B["M_GG"]
    E={tau:expm(tau*M) for tau in (0.02,0.25,1.0)}
    for gam in (0.2,0.3,0.4,0.45,0.5):
        w=np.exp(-gam*z)
        out=[]
        for tau,Et in E.items():
            A=(w[:,None]*Et)/w[None,:]
            out.append(f"tau={tau}: {np.linalg.norm(A,2)*np.exp(2*gam*tau):.4f}")
        print(f"N={N} zmin={zmin} gamma={gam}: ||e^(tau M)||_gamma * e^(2 gamma tau): "+"  ".join(out)); sys.stdout.flush()
