import numpy as np, sys, time
from common import *
from mycolloc import Qmat, csym
gams=(0.1,0.2,0.25,0.26,0.27,0.3,0.4,0.5)
cases=eval(sys.argv[1]); zmax=float(sys.argv[2]) if len(sys.argv)>2 else 7.0
print("zmax",zmax,"cols: gamma | repo new 0.5*Q_op | repo new M_GG+2Dz | repo old-diag(tail) | my collocation | ratio repo/mine")
for (N,zmin) in cases:
    t0=time.time()
    Bn=NEW.build_blocks_z(N=N,zmin=zmin,zmax=zmax,edge_bc=True)
    Bo=NEW.build_blocks_z(N=N,zmin=zmin,zmax=zmax,edge_bc=True,hilbert_grid_support=False)
    z=Bn['z']; Fq=0.5*Bn['Q_op']; Fm=Bn['M_GG']+2*Bn['Dz']; Fo=0.5*Bo['Q_op']
    zf=np.linspace(zmin,zmax,N); Qm=Qmat(zf)[1:,1:]; Fc=0.5*Qm
    print(f"N={N} zmin={zmin}  (build {time.time()-t0:.0f}s)  ||Q_repo-Q_mine||_2/||Q_mine||_2 = {np.linalg.norm(Bn['Q_op']-Qm,2)/np.linalg.norm(Qm,2):.3e}; max|forcing-0.5Q_op|={abs(Fq-Fm).max():.2e}",flush=True)
    for g in gams:
        a,b,o,c=csym(Fq,z,g),csym(Fm,z,g),csym(Fo,z,g),csym(Fc,z,g)
        print(f"   {g:4.2f} | {a:+.6e} | {b:+.6e} | {o:+.6e} | {c:+.6e} | {a/c:8.4f}",flush=True)
