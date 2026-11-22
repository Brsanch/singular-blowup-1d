"""Domain / grid / sponge convergence of the reviewer's solver, and how clean the gate domain [-64,64] n=4096 is."""
import numpy as np, sys
from lib import *
from solver import *
from scipy.interpolate import CubicSpline
sc = G0_maxnorm_scale()
A=float(sys.argv[1]) if len(sys.argv)>1 else 0.05
def run(zmin,zmax,n,nst,sponge=None,S=1.5):
    s=Solver(zmin,zmax,n,sponge=sponge); h=A*h0_family(s)/sc
    H,bad=s.run(h,S,nst,nst//15)
    return s,H,bad
s0,H0,_=run(-128,128,16384,600,sponge=16.0)
G0=np.array([s0.N(h) for h in H0]); zz=s0.g.z
c1=[np.mean((G0[j]/s0.g.w)[(zz>-18)&(zz<-15)]) for j in range(len(G0))]
print("reference: [-128,128] n=16384 sponge 16; c1(s)=G/w at z in (-18,-15):", ["%.4e"%c for c in c1])
print("config: max over region of |G - G_ref| at s = 0.1, 0.5, 1.0, 1.5  (edge column is |G/w - ref|)")
for (zmin,zmax,n,sp) in ((-64,64,4096,None),(-48,48,4096,None),(-36,36,4096,None),(-72,72,4096,None),(-72,72,8192,None),(-96,96,8192,None),(-96,96,16384,None),(-104,24,4096,None),(-64,64,4096,10.0),(-64,64,8192,10.0),(-96,96,16384,12.0),(-104,24,8192,8.0)):
    s,H,bad=run(zmin,zmax,n,600,sponge=sp)
    tag=f"[{zmin},{zmax}] n={n} sponge={sp}"
    nj=len(H)
    out=[]
    for j in (1,5,10,15):
        if j>=nj: out.append("NaN(s=%.3f)"%bad); continue
        G=s.N(H[j]); cs=CubicSpline(s.g.z,G); row=[]
        for lo,hi in ((-18,-12),(-12,12),(12,24)):
            mk=(zz>lo)&(zz<hi)&(zz>zmin+1)&(zz<zmax-1)
            e=np.max(np.abs(cs(zz[mk])-G0[j][mk])/(np.exp(zz[mk]) if hi<=-12 else 1)); row.append("%.1e"%e)
        out.append("/".join(row))
    print(f"  {tag:36s} edge/bulk/far: "+"   ".join(out))
