import numpy as np
from common import *; from exactmat import exact_Q
for N,zmin in ((400,-10.),(800,-10.),(1600,-10.),(800,-4.),(800,-14.)):
    z=np.linspace(zmin,7.,N); Qf=exact_Q(z); Bf=S.build_blocks_z(N=N,zmin=zmin,zmax=7.,edge_bc=False)
    M=-2*Bf['Dz']+0.5*Qf
    v,V=np.linalg.eig(M); j=np.argmax(v.real); e=np.abs(V[:,j]); 
    nplus=v[v.real>0]
    print(f"N={N} zmin={zmin}: exact-kernel free-edge rightmost {v[j]:.4f}; all Re>0: {np.round(nplus,3)}; eigvec |V| max at z={z[e.argmax()]:.2f}, mass within 0.5 of zmin={np.sum(e[z<zmin+0.5]**2)/np.sum(e**2):.3f}; Q[0,0]/2={0.5*Qf[0,0]:+.3f} 3/dz={3/(z[1]-z[0]):.1f}; repo-tail: Q[0,0]/2={0.5*Bf['Q_op'][0,0]:+.3f}")
    # same but hat-Galerkin-free alternative closure: upwind (first-order) transport to test dependence on Dz closure
    n=N; dz=z[1]-z[0]; Du=np.zeros((n,n))
    for i in range(1,n): Du[i,i]=1/dz; Du[i,i-1]=-1/dz
    Du[0,0]=-1/dz; Du[0,1]=1/dz
    for nm,Q in (('exact',Qf),('repo-tail',Bf['Q_op'])):
        vv=np.linalg.eigvals(-2*Du+0.5*Q); print(f"      upwind-transport closure, {nm}: maxRe={vv.real.max():+.4f} n+={int((vv.real>0).sum())}")
