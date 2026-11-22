# C*(sigma) := ||S(tau)|| e^{2 sigma tau}  (tau-independent by scaling) -- lower bounds from enlarging low-frequency families
import numpy as np, sys
from a4 import run
taus=[0.05,0.7,3.0]
for sig in [float(x) for x in sys.argv[1:]]:
    row=[]
    for L in [4,8,16,32]:
        dt=min(0.02,np.pi/(2*L+10)/4)
        n,k,out=run(sig,[1.0,0.2],np.exp(np.arange(-L,L+0.01,0.5 if L<=16 else 1.0)),[0,0.35,-0.35],taus,T=40.,dt=dt,thr=1e-10)
        row.append((L,n,k,[o[2] for o in out]))
    print(f"sigma={sig}: "+" | ".join(f"L={L} n={n} k={k}: "+",".join(f"{v:.6f}" for v in vs) for L,n,k,vs in row)); sys.stdout.flush()
