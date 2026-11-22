import numpy as np, sys
from lib import *
from solver import *
from scipy.interpolate import CubicSpline
sc = G0_maxnorm_scale()
def run(zmin,zmax,n,nst,A,sponge=None,filt=None,S=1.5):
    s=Solver(zmin,zmax,n,filt=filt,sponge=sponge); h=A*h0_family(s)/sc
    H,bad=s.run(h,S,nst,nst//10)
    return s,H,bad
ref=None
for (zmin,zmax,n,sp) in ((-64,64,4096,None),(-64,64,4096,10.0),(-64,64,8192,10.0),(-64,64,16384,10.0),(-96,96,8192,12.0),(-96,96,16384,12.0),(-96,96,16384,None),(-48,48,4096,8.0),(-104,24,8192,8.0)):
    s,H,bad=run(zmin,zmax,n,600,0.05,sponge=sp)
    if bad: print(zmin,zmax,n,sp,"NaN at",bad); continue
    G=np.array([s.N(h) for h in H])
    if ref is None: ref=(s,G); print(zmin,zmax,n,sp,"reference"); continue
    zz=ref[0].g.z; out=[]
    for lo,hi in ((-20,-12),(-12,12),(12,24)):
        mk=(zz>lo)&(zz<hi); e=0
        for j in range(len(G)):
            e=max(e,np.max(np.abs(CubicSpline(s.g.z,G[j])(zz[mk])-ref[1][j][mk])/(np.exp(zz[mk]) if hi<=-12 else 1)))
        out.append("%.2e"%e)
    print(zmin,zmax,n,sp,"finite; diff vs ref: edge(G/w) (-20,-12), bulk (-12,12), far (12,24):",out)
