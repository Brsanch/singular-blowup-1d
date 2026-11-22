import numpy as np
from common import *
from mycolloc import csym
orig=NEW.hilbert_support_z
# (a) with edge_bc=True the result does not depend on the row-0 diagonal
def Hmod(Y,Ymax,grid_support=False):
    M=orig(Y,Ymax,grid_support=grid_support); M[0,0]+=123.456; return M
B1=NEW.build_blocks_z(N=500,zmin=-12,zmax=7,edge_bc=True)
NEW.hilbert_support_z=Hmod; B2=NEW.build_blocks_z(N=500,zmin=-12,zmax=7,edge_bc=True); NEW.hilbert_support_z=orig
print("(a) edge_bc=True: all returned arrays identical when Hz[0,0] is perturbed by +123:",all(np.array_equal(np.asarray(B1[k]),np.asarray(B2[k])) for k in B1))
# (b) gate (1) closed form H((Y-1)^{-1/2})=0, historical vs grid-supported diagonal
for N,zmin,zmax in ((1600,-14.,7.),(3200,-14.,7.),(1600,-6.,7.)):
    z=np.linspace(zmin,zmax,N);Y=1+np.exp(z);Ymax=Y[-1]+.5*(Y[-1]-Y[-2]);g=(Y-1)**-.5;win=(Y>1.05)&(Y<300)
    for gs in (False,True):
        Hg=orig(Y,Ymax,grid_support=gs)@g
        print(f"(b) N={N} zmin={zmin} grid_support={gs}: max|H(g)| on Y in (1.05,300) = {abs(Hg[win]).max():.3e};  max |H(g)|*sqrt(v) over all nodes>=1: {abs(Hg[1:]*np.sqrt(Y[1:]-1)).max():.3e}")
# (c) top-end mismatch: rows near the top, bump at zmax-1; exact for g supported on the grid
from exact import exact_Hg,bump
for N in (800,3200):
    zmin,zmax=-10.,7.;z=np.linspace(zmin,zmax,N);Y=1+np.exp(z);Ymax=Y[-1]+.5*(Y[-1]-Y[-2])
    G=lambda t: np.exp(-(zmax-np.asarray(t)));Gv=G(z)   # nonzero at the top node (jump at z_max)
    Hn=orig(Y,Ymax,grid_support=True)@(np.exp(-z/2)*Gv)
    idx=[N-2,N-3,N-5,N-10,int(np.argmin(abs(z-(zmax-.5)))),int(np.argmin(abs(z-(zmax-1)))),int(np.argmin(abs(z-(zmax-3))))]
    ex=np.array([exact_Hg(G,z[i],zmin,zmax) for i in idx])
    print(f"(c) N={N} top rows z_max-z={np.round(zmax-z[idx],3)} rel err Hg (G=e^(z-zmax), jump at z_max):",np.round((Hn[idx]-ex)/abs(ex),4))
# (d) floor at gamma<=1/4: N dependence incl. 400 and 6400, and with the top node removed
for N in (400,800,1600,3200,6400):
    B=NEW.build_blocks_z(N=N,zmin=-10.,zmax=12.,edge_bc=True);F=.5*B['Q_op'];z=B['z']
    print(f"(d) zmin=-10 zmax=12 N={N}: C(0.25)={csym(F,z,.25):+.4e}  C(0.10)={csym(F,z,.10):+.4e}  top 3 nodes removed: {csym(F[:-3,:-3],z[:-3],.25):+.4e}; top 10% removed: {csym(F[:-(N//10),:-(N//10)],z[:-(N//10)],.25):+.4e}",flush=True)
