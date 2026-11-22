import numpy as np, sys
from a4 import run
taus=[1e-3,1e-2,0.03,0.1,0.2,0.4,0.7,1.0,1.5,2.5,4.0]
sig=float(sys.argv[1])
fam=sys.argv[2]
F={"medium":([1.0],np.exp(np.linspace(-6,4,11)),[0,0.6,-0.6,1.2,-1.2,1.45,-1.45]),
   "wide":([1.0],np.exp(np.linspace(-12,6,37)),[0,0.7,-0.7,1.3,-1.3]),
   "lowfreq":([1.0,0.2],np.exp(np.linspace(-16,8,49)),[0,0.3,-0.3]),
   "hf":([1.0],np.exp(np.linspace(-6,4,11)),[0,1.0,-1.0,1.4,-1.4,1.5,-1.5,1.54,-1.54,1.56,-1.56])}
a,r,b=F[fam]
T,dt=(300.,0.04) if fam!="hf" else (1500.,0.05)
n,k,out=run(sig,a,r,b,taus,T=T,dt=dt)
print(f"sigma={sig} family={fam} n={n} rank_kept={k}")
print("   tau: ratio | ratio*e^{2 sig tau} | best single member: "+"; ".join(f"{o[0]:g}: {o[1]:.5f} | {o[2]:.5f} | {o[3]:.4f}" for o in out))
