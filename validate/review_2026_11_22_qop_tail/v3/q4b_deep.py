import numpy as np, sys
from common import *
from mycolloc import Qmat, csym
gams=(0.1,0.2,0.25,0.252,0.255,0.26,0.27,0.3,0.4)
for (N,zmin,zmax) in eval(sys.argv[1]):
    Bn=NEW.build_blocks_z(N=N,zmin=zmin,zmax=zmax,edge_bc=True); z=Bn['z']; Fq=0.5*Bn['Q_op']
    Fc=0.5*Qmat(np.linspace(zmin,zmax,N))[1:,1:]
    print(f"N={N} zmin={zmin} zmax={zmax} dz={Bn['dz']:.4f}")
    for g in gams:
        a,c=csym(Fq,z,g),csym(Fc,z,g)
        # where does the top eigenvector live?
        sw=np.exp(-g*z); At=(sw[:,None]*Fq)/sw[None,:]; w,V=np.linalg.eigh(0.5*(At+At.T)); v=V[:,-1]; zc=(z*v**2).sum()
        print(f"   {g:5.3f} | repo {a:+.6e} | mine {c:+.6e} | ratio {a/c:8.4f} | repo top eigvec centroid z={zc:+.2f}",flush=True)
