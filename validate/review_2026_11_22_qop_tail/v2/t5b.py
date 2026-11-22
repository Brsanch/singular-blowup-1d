import numpy as np, sys
from multiprocessing import Pool
from common import *; from exactmat import exact_Q
def lam(F,z,g,mask=None):
    w=np.exp(-g*z); A=(w[:,None]*F)/w[None,:]; Sm=0.5*(A+A.T)
    if mask is not None: Sm=Sm[np.ix_(mask,mask)]
    return np.linalg.eigvalsh(Sm)
def ramp_A(Q,z,dz,gamma=0.4,z0=-6.,s=2.):
    pp=2*gamma*0.5*(1+np.tanh((z-z0)/s)); phi=np.cumsum(pp)*dz; sw=np.exp(-0.5*phi)
    Ft=(sw[:,None]*(0.5*Q))/sw[None,:]
    return np.linalg.eigvalsh(-np.diag(pp)+0.5*(Ft+Ft.T))[-1]
def run(args):
    N,zmin,zmax=args; out=[]
    z=np.linspace(zmin,zmax,N); Qx=exact_Q(z)[1:,1:]
    B=S.build_blocks_z(N=N,zmin=zmin,zmax=zmax,edge_bc=True); zz=B['z']; Fr=0.5*B['Q_op']; Fx=0.5*Qx
    s=f"N={N} zmin={zmin} zmax={zmax} bc=1\n"
    for g in (0.20,0.24,0.25,0.26,0.27,0.28,0.29,0.30,0.32,0.40):
        er=lam(Fr,zz,g); ex=lam(Fx,zz,g)
        m=zz<zmin+8
        er2=lam(Fr,zz,g,m); ex2=lam(Fx,zz,g,m)
        s+=f"  g={g:.2f} repo lmax={er[-1]:+.3e} (lmin={er[0]:+.2e})  exact lmax={ex[-1]:+.3e} (lmin={ex[0]:+.2e}) | edge block z<zmin+8: repo {er2[-1]:+.3e} exact {ex2[-1]:+.3e}\n"
    s+=f"  ramp A_full(gamma=.4,z0=-6,s=2): repo {ramp_A(B['Q_op'],zz,B['dz']):+.4f}  exact-kernel {ramp_A(Qx,zz,B['dz']):+.4f}\n"
    return s
if __name__=="__main__":
    cfg=[(400,-10.,12.),(800,-10.,12.),(1600,-10.,12.),(3200,-10.,12.),(1400,-12.,12.),(1400,-14.,12.),(1400,-16.,12.),(900,-12.,12.),(2800,-12.,12.)]
    with Pool(4) as p:
        for o in p.imap(run,cfg): print(o); sys.stdout.flush()
