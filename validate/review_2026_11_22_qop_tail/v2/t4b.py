import numpy as np, sys
from multiprocessing import Pool
from common import *; from exactmat import exact_Q
def one(a):
    N,zmin,zmax=a
    B=S.build_blocks_z(N=N,zmin=zmin,zmax=zmax,edge_bc=False,edge_tail=True)
    v=np.linalg.eigvals(B['M_GG'])
    rl=v[np.abs(v.imag)<1e-7].real
    j=np.argmax(v.real)
    return f"N={N:4d} zmin={zmin:5.1f} zmax={zmax:4.1f} free+tail: maxRe={v[j].real:+.6f} (Im {v[j].imag:+.1f}) n(Re>0)={int((v.real>0).sum())} rightmost REAL eig={rl.max() if len(rl) else np.nan:+.5f} low(|Im|<20)={v[np.abs(v.imag)<20].real.max():+.4f}"
def stiff(a):
    N,zmin,zmax=a
    z=np.linspace(zmin,zmax,N); Qx=exact_Q(z)[1:,1:]
    B=S.build_blocks_z(N=N,zmin=zmin,zmax=zmax,edge_bc=True); B0=S.build_blocks_z(N=N,zmin=zmin,zmax=zmax,edge_bc=True,edge_tail=False)
    s=f"N={N} zmin={zmin} zmax={zmax} bc=1:"
    for nm,M in (('repo tail',B['M_GG']),('truncated',B0['M_GG']),('exact-kernel',-2*B['Dz']+0.5*Qx)):
        v=np.linalg.eigvals(M); j=np.argmin(v.real); 
        s+=f"\n    {nm:12s} minRe eig={v[j].real:+.2f}{v[j].imag:+.2e}i  maxRe={v.real.max():+.5f} n+={int((v.real>0).sum())} low(|Im|<20)={v[np.abs(v.imag)<20].real.max():+.4f} rho={np.abs(v).max():.1f}"
    Qf=exact_Q(z); Bf=S.build_blocks_z(N=N,zmin=zmin,zmax=zmax,edge_bc=False)
    v=np.linalg.eigvals(-2*Bf['Dz']+0.5*Qf); 
    s+=f"\n    FREE-EDGE exact-kernel: maxRe={v.real.max():+.5f} n+={int((v.real>0).sum())} low(|Im|<20)={v[np.abs(v.imag)<20].real.max():+.4f} minRe={v.real.min():+.1f}"
    return s
if __name__=="__main__":
    cfg=[(N,float(z),zm) for z in (0,-1,-2,-3) for (N,zm) in ((200,7.),(800,7.),(800,16.),(1600,24.))]+[(50,-6.,7.),(50,-14.,7.),(30,-10.,7.),(1600,-4.,24.),(2400,-6.,30.),(2000,-20.,7.),(2000,-24.,7.)]
    with Pool(4) as p:
        for o in p.imap(one,cfg): print(o); sys.stdout.flush()
        for o in p.imap(stiff,[(800,-10.,7.),(800,-14.,7.),(1600,-14.,7.),(800,-6.,7.)]): print(o); sys.stdout.flush()
