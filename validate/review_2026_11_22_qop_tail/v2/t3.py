import numpy as np, sys
from multiprocessing import Pool
from common import *
def one(a):
    N,zmin,zmax,bc=a; out=[]
    r={}
    for tail in (False,True):
        B=S.build_blocks_z(N=N,zmin=zmin,zmax=zmax,edge_bc=bc,edge_tail=tail)
        v=np.linalg.eigvals(B['M_GG']); r[tail]=v
        dt=0.25*B['dz']
        rho=np.abs(v).max()
        realv=v[np.abs(v.imag)<1e-8].real
        low=v[np.abs(v.imag)<20]
        out.append(f"N={N:4d} zmin={zmin:5.1f} zmax={zmax:4.1f} bc={int(bc)} tail={int(tail)}: maxRe={v.real.max():+.5f} #Re>1e-6={int((v.real>1e-6).sum())} maxRe(|Im|<20)={low.real.max():+.4f} most-neg-real={realv.min() if len(realv) else np.nan:+9.2f} minRe={v.real.min():+9.2f} rho={rho:8.2f} dt*rho={dt*rho:.3f} dt*|minRe|={dt*abs(v.real.min()):.3f}  -0.13e^(-zmin/2)={-0.13*np.exp(-zmin/2):+.1f}")
    return "\n".join(out)
if __name__=="__main__":
    cfg=[(800,-10.,7.,True),(1600,-10.,7.,True),(800,-14.,7.,True),(1600,-14.,7.,True),(800,-16.,7.,True),(1600,-16.,7.,True),
         # marching grids
         (1400,-14.,12.,True),(1100,-12.,12.,True),(1800,-16.,14.,True),(1400,-12.,12.,True),(1000,-10.,12.,True),(1800,-14.,12.,True),(500,-12.,8.,True),(1100,-12.,10.,True),
         (1000,-10.,12.,False),(1400,-12.,12.,False),(1800,-14.,12.,False),(2200,-16.,12.,False),
         (1400,-14.,7.,True),(1000,-14.,8.,True),(1000,-14.,10.,True)]
    with Pool(4) as p:
        for o in p.imap(one,cfg): print(o); sys.stdout.flush()
