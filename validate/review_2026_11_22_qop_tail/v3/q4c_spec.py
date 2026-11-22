import numpy as np, sys
from common import *
from mycolloc import Qmat
def summ(M):
    v=np.linalg.eigvals(M); k=np.argmin(v.real)
    lo=v[abs(v.imag)<20]
    return f"maxRe={v.real.max():+.5f} #Re>0={int((v.real>0).sum())} #Re>1e-6={int((v.real>1e-6).sum())} right|Im|<20={lo.real.max():+.4f} leftmost={v[k].real:+.2f}{v[k].imag:+.2f}i rho={abs(v).max():.1f}"
for (N,zmin,zmax) in eval(sys.argv[1]):
    print(f"N={N} zmin={zmin} zmax={zmax}")
    for tag,kw in (("truncated      ",dict(edge_tail=False)),("tail+hist (ad4)",dict(edge_tail=True,hilbert_grid_support=False)),("default new    ",dict())):
        B=NEW.build_blocks_z(N=N,zmin=zmin,zmax=zmax,edge_bc=True,**kw)
        print("   ",tag,summ(B['M_GG']),flush=True)
    Qm=Qmat(np.linspace(zmin,zmax,N))[1:,1:]
    print("    my colloc + Dz ",summ(-2*B['Dz']+0.5*Qm),flush=True)
