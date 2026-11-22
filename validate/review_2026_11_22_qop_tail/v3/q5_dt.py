import numpy as np
from common import *
R=lambda x:1+x+x**2/2+x**3/6+x**4/24
grids=[("nonlinear_shape_march main",1400,-14.,12.,True),("nsm conv",1100,-12.,12.,True),("nsm conv",1800,-16.,14.,True),
       ("coevolve main (mf/stress/proj/paper_match)",1400,-12.,12.,False),("coevolve conv",1000,-10.,12.,False),("coevolve conv",1800,-14.,12.,False),
       ("coevolve_projected step0",300,-12.,12.,False),("coevolve_projected step0",600,-12.,12.,False),("coevolve_projected step0",1000,-12.,12.,False),
       ("fig_coevolve",500,-12.,8.,False),("fig_nonlinear (uses nsm.march)",1100,-12.,10.,True)]
print("grid | variant | dt*rho at dt=0.25dz | dt used | dt_used*rho | max|R(dt lam)| full M_GG | same with first node pinned | worst eig")
for name,N,zmin,zmax,cap in grids:
    for tag,kw in (("ad4 default(tail+hist)",dict(hilbert_grid_support=False)),("NEW default",dict())):
        B=NEW.build_blocks_z(N=N,zmin=zmin,zmax=zmax,edge_bc=True,**kw); dz=B['dz']
        v=np.linalg.eigvals(B['M_GG']); v2=np.linalg.eigvals(B['M_GG'][1:,1:]); rho=abs(v).max()
        dt0=0.25*dz; dt=min(dt0,2/rho) if cap else dt0
        g=abs(R(dt*v)); g2=abs(R(dt*v2)); k=np.argmax(g)
        print(f"{name:42s} N={N} zmin={zmin} zmax={zmax} | {tag:22s} | {dt0*rho:6.3f} | dt={dt:.5f} | {dt*rho:6.3f} | {g.max():.6f} | {g2.max():.6f} | {v[k]:.1f}  (hh block: {abs(R(dt*np.linalg.eigvals(B['M_hh']))).max():.6f})",flush=True)
