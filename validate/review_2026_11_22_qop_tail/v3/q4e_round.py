import numpy as np
from common import *
from mycolloc import Qmat, csym
N=1600;zmax=12.
for zmin in (-22.,-24.,-26.,-28.,-30.):
    z=np.linspace(zmin,zmax,N);Y=1+np.exp(z);v=np.exp(z)
    rel=abs((Y[1:]-Y[:-1])-(v[1:]-v[:-1]))/(v[1:]-v[:-1])
    B=NEW.build_blocks_z(N=N,zmin=zmin,zmax=zmax,edge_bc=True);F=.5*B['Q_op'];zz=B['z']
    Fc=.5*Qmat(z)[1:,1:]
    print(f"zmin={zmin}: max rel error of node spacing from Y=1+e^z: {rel.max():.1e} | C(0.25): repo {csym(F,zz,.25):+.3e} (top3 removed {csym(F[:-3,:-3],zz[:-3],.25):+.3e}) mine {csym(Fc,zz,.25):+.3e} | C(0.10): repo {csym(F,zz,.10):+.3e} mine {csym(Fc,zz,.10):+.3e}",flush=True)
