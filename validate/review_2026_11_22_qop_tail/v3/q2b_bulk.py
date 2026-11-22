import numpy as np
from common import *
from exact import *
zmax=7.0
F={'G=1/(1+v)  (G(edge)=1)':lambda zp: 1/(1+np.exp(zp)),
   'G=v^.5/(1+v)^1.5 (g~const at edge)':lambda zp: np.exp(.5*zp)/(1+np.exp(zp))**1.5,
   'G=v^.25/(1+v)':lambda zp: np.exp(.25*zp)/(1+np.exp(zp))}
for name,G in F.items():
  for zmin in (-6.,-10.,-14.):
    for N in (800,3200):
        Bo=NEW.build_blocks_z(N=N,zmin=zmin,zmax=zmax,edge_bc=False,edge_tail=True,hilbert_grid_support=False)
        Bn=NEW.build_blocks_z(N=N,zmin=zmin,zmax=zmax,edge_bc=False,edge_tail=True)
        Bt=NEW.build_blocks_z(N=N,zmin=zmin,zmax=zmax,edge_bc=False,edge_tail=False)
        z=Bo['z']; idx=sorted(set([1,2,8]+[int(np.argmin(abs(z-t))) for t in (zmin+.5,zmin+1,zmin+2,zmin+3,-2.,0.,2.)]))
        Gv=G(z)
        exQ=np.array([exact_Q(G,z[i],-90,90) for i in idx]); exH=np.array([exact_Hg(G,z[i],-90,90) for i in idx])
        print(f"\n{name} zmin={zmin} N={N}\n   z-zmin: "+" ".join(f"{z[i]-zmin:9.3f}" for i in idx))
        print("   exact Q "+" ".join(f"{x:9.5f}" for x in exQ))
        for tag,B in (('trunc/hist',Bt),('tail/hist ',Bo),('tail/grid ',Bn)):
            print(f"   relerr Q  {tag} "+" ".join(f"{x:9.2e}" for x in ((B['Q_op']@Gv)[idx]-exQ)/abs(exQ)))
        for tag,B in (('hist',Bo),('grid',Bn)):
            print(f"   relerr Hg {tag} "+" ".join(f"{x:9.2e}" for x in ((B['Hg_op']@Gv)[idx]-exH)/abs(exH)))
