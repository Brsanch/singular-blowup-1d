"""Look at |Ghat(1/4 + i xi, s)| of the nonlinear solution under different numerical settings."""
import numpy as np, sys, math
from lib import *
from solver import *
sc = G0_maxnorm_scale()
XI=(2,4,6,8,10,12,14,16,18,20,24,28,32,40)
def spec(zmin,zmax,n,nst,A,sponge,sstop=0.5,filt=None):
    s=Solver(zmin,zmax,n,sponge=sponge,filt=filt); h=A*h0_family(s)/sc
    H,bad=s.run(h,sstop,nst,nst)
    if bad: return None
    G=s.N(H[-1]); Gh=np.abs(s.g.hat(G,0.25))/A
    hh=np.abs(s.g.hat(H[-1],0.25))/A
    return [Gh[int(round(x/s.g.dxi))] for x in XI],[hh[int(round(x/s.g.dxi))] for x in XI]
print("xi:           "+" ".join(f"{x:9d}" for x in XI))
for lab,args in (("A=.05 [-64,64] 8192 sp10 ds=.5/400",(-64,64,8192,400,0.05,10.0)),
                 ("A=.05 [-64,64] 8192 sp10 ds=.5/1600",(-64,64,8192,1600,0.05,10.0)),
                 ("A=.05 [-64,64] 16384 sp10",(-64,64,16384,400,0.05,10.0)),
                 ("A=.05 [-64,64] 8192 sp20",(-64,64,8192,400,0.05,20.0)),
                 ("A=.05 [-128,128] 16384 none",(-128,128,16384,400,0.05,None)),
                 ("A=.05 [-104,24] 8192 sp8",(-104,24,8192,400,0.05,8.0)),
                 ("A=.005 [-64,64] 8192 sp10",(-64,64,8192,400,0.005,10.0)),
                 ("A=.0005 [-64,64] 8192 sp10",(-64,64,8192,400,0.0005,10.0)),
                 ("A=1e-9 (linear) ",(-64,64,8192,400,1e-9,10.0)),
                 ("A=.2 [-64,64] 8192 sp10",(-64,64,8192,400,0.2,10.0)),
                 ):
    r=spec(*args)
    if r is None: print(lab,"NaN"); continue
    print(f"{lab:36s} |Ghat|/A: "+" ".join(f"{v:9.2e}" for v in r[0]))
    print(f"{'':36s} |hhat|/A: "+" ".join(f"{v:9.2e}" for v in r[1]))
