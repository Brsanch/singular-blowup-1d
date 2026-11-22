import numpy as np
from common import *; from fastq import *
zmin=-14.;G=Ga
for z in (-14.,-13.9,-10.,0.,5.):
    print(z, Q_exact(G,z,zmin,7.), Q_fast(G,[z],zmin,7.)[0], Q_fast(G,[z],-60,40.)[0])
G=mk_bump(zmin+1)
for z in (-14.,-13.2,-13.,-10.,0.):
    print(z, Q_exact(G,z,zmin,7.,(zmin+1,)), Q_fast(G,[z],zmin,7.)[0])
print("a exact", -1/(4*np.sqrt(2)))
