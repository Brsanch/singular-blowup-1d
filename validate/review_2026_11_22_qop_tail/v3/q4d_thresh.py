import numpy as np
from common import *
from mycolloc import Qmat, csym
gs=np.round(np.arange(0.250,0.2721,0.002),3)
for (N,zmin,zmax) in [(1600,-10.,12.),(1600,-10.,7.),(1600,-14.,12.),(1600,-22.,12.),(1600,-30.,12.)]:
    B=NEW.build_blocks_z(N=N,zmin=zmin,zmax=zmax,edge_bc=True);F=.5*B['Q_op'][:-3,:-3];z=B['z'][:-3]
    Fc=.5*Qmat(np.linspace(zmin,zmax,N))[1:,1:]; zc=B['z']
    r=[csym(F,z,g) for g in gs]; c=[csym(Fc,zc,g) for g in gs]
    fp=lambda a: next((g for g,x in zip(gs,a) if x>0),None)
    print(f"N={N} zmin={zmin} zmax={zmax}: first gamma with value>0: repo(top 3 nodes removed) {fp(r)}, my collocation {fp(c)}")
    print("    repo-3:",[f"{x:+.1e}" for x in r]); print("    mine  :",[f"{x:+.1e}" for x in c],flush=True)
