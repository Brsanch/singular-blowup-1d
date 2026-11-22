import numpy as np, sys
from common import *
from mycolloc import Qmat
orig=NEW.hilbert_support_z
def make(variant):
    def H(Y,Ymax,grid_support=False):
        M=orig(Y,Ymax,grid_support=grid_support)
        if grid_support and variant!="hist":
            n=len(Y); wY0=0.5*(Y[1]-Y[0]); K0=(np.r_[0.5*(Y[1]-Y[0]),0.5*(Y[2:]-Y[:-2]),0.5*(Y[-1]-Y[-2])]/(Y[0]-Y+np.eye(1,n,0)[0]))[1:].sum()
            if variant=="halfcell": pv0=np.log(wY0/(Ymax-Y[0]))          # g=g_0 on half a cell below Y_0 (mirror of the top end)
            elif variant=="cellavg": pv0=np.log(wY0)-1-np.log(Ymax-Y[0]) # average of the exact pv over [Y_0,Y_0+wY0]
            M[0,0]=(pv0-K0)/np.pi
        return M
    return H
def top(M,k=3):
    v=np.linalg.eigvals(M); i=np.argsort(-v.real)[:k]; return "  ".join(f"{v[j].real:+.4f}{v[j].imag:+.2f}i" for j in i)
zmax=7.0
for zmin in (-12.,-6.):
  for N in (400,800,1000,1600,3200):
    print(f"zmin={zmin} N={N}")
    NEW.hilbert_support_z=orig
    Bt=NEW.build_blocks_z(N=N,zmin=zmin,zmax=zmax,edge_bc=False,edge_tail=False); print("   truncated/hist       :",top(Bt['M_GG']))
    Bh=NEW.build_blocks_z(N=N,zmin=zmin,zmax=zmax,edge_bc=False,edge_tail=True,hilbert_grid_support=False); print("   tail/hist (ad4)      :",top(Bh['M_GG']))
    for var in ("hist","halfcell","cellavg"):
        NEW.hilbert_support_z=make(var)
        B=NEW.build_blocks_z(N=N,zmin=zmin,zmax=zmax,edge_bc=False,edge_tail=True)
        v,V=np.linalg.eig(B['M_GG']); k=np.argmax(v.real); e=abs(V[:,k]); 
        print(f"   tail/grid row0={var:8s}:",top(B['M_GG']),f"| Hg[0,0]*sqrt(v0)*pi={B['Hg_op'][0,0]*np.exp(zmin/2)*np.pi:+.3f} forcing[0,0]={(B['M_GG']+2*B['Dz'])[0,0]:+.3f} | top mode |V| at nodes 0..3 / max: {np.round(e[:4]/e.max(),3)}",flush=True)
    NEW.hilbert_support_z=orig
    Qm=Qmat(np.linspace(zmin,zmax,N)); print("   my exact collocation :",top(-2*Bt['Dz']+0.5*Qm), f"| 0.5*Q[0,0]={0.5*Qm[0,0]:+.3f}",flush=True)
