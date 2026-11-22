import numpy as np, inspect
from common import *
print("PRE sig:",inspect.signature(PRE.build_blocks_z)); print("AD4 sig:",inspect.signature(AD4.build_blocks_z)); print("NEW sig:",inspect.signature(NEW.build_blocks_z))
bad=0;n=0
for (N,zmin,zmax) in [(60,-6,7),(300,-10,7),(500,-12,8),(700,-14,12),(400,-18,7),(333,-3,5)]:
  for bc in (True,False):
    for beta in (0.0,-0.25,-0.5,0.3):
      a=NEW.build_blocks_z(N=N,zmin=zmin,zmax=zmax,beta=beta,edge_bc=bc,edge_tail=False)
      b=AD4.build_blocks_z(N=N,zmin=zmin,zmax=zmax,beta=beta,edge_bc=bc,edge_tail=False)
      kw={} if 'edge_tail' not in inspect.signature(PRE.build_blocks_z).parameters else dict(edge_tail=False)
      c=PRE.build_blocks_z(N=N,zmin=zmin,zmax=zmax,beta=beta,edge_bc=bc,**kw)
      for k in a:
        n+=1
        if k not in b or not np.array_equal(np.asarray(a[k]),np.asarray(b[k])): bad+=1; print("DIFF vs AD4",N,zmin,bc,beta,k)
        if k not in c or not np.array_equal(np.asarray(a[k]),np.asarray(c[k])): bad+=1; print("DIFF vs PRE",N,zmin,bc,beta,k)
      # explicit hilbert_grid_support=False with tail on == AD4 tail on
      a2=NEW.build_blocks_z(N=N,zmin=zmin,zmax=zmax,beta=beta,edge_bc=bc,edge_tail=True,hilbert_grid_support=False)
      b2=AD4.build_blocks_z(N=N,zmin=zmin,zmax=zmax,beta=beta,edge_bc=bc,edge_tail=True)
      for k in a2:
        n+=1
        if not np.array_equal(np.asarray(a2[k]),np.asarray(b2[k])): bad+=1; print("DIFF tail-on/hist vs AD4 tail-on",N,zmin,bc,beta,k)
print("comparisons",n,"mismatches",bad)
# direct hilbert_support_z default
z=np.linspace(-10,7,300);Y=1+np.exp(z);Ymax=Y[-1]+.5*(Y[-1]-Y[-2])
print("hilbert default identical to AD4/PRE:",np.array_equal(NEW.hilbert_support_z(Y,Ymax),AD4.hilbert_support_z(Y,Ymax)),np.array_equal(NEW.hilbert_support_z(Y,Ymax),PRE.hilbert_support_z(Y,Ymax)))
Hn=NEW.hilbert_support_z(Y,Ymax,grid_support=True);Ho=NEW.hilbert_support_z(Y,Ymax)
D=Hn-Ho; off=D-np.diag(np.diag(D)); print("offdiag change max",abs(off).max(),"row0 diag change",D[0,0])
ex=np.log((Y-Y[0])[1:]/(Ymax-Y)[1:])/np.pi
print("H@1 rows>=1 vs exact H[1_[Y0,Ymax]]: max abs err",abs((Hn@np.ones(len(Y)))[1:]-ex).max())
print("row0 H@1 new:",(Hn@np.ones(len(Y)))[0],"= ln((Y0-1)/(Ymax-Y0))/pi =",np.log((Y[0]-1)/(Ymax-Y[0]))/np.pi)
