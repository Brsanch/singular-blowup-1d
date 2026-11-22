import numpy as np, sys
from multiprocessing import Pool
from common import *; from exactmat import exact_Q
gam=np.arange(0.05,0.5001,0.05)
def lam(F,z,g):
    w=np.exp(-g*z); A=(w[:,None]*F)/w[None,:]
    ev,V=np.linalg.eigh(0.5*(A+A.T)); return ev[-1],V[:,-1]
def run(args):
    N,zmin,zmax=args; out=[]
    z=np.linspace(zmin,zmax,N); Qx=exact_Q(z)
    for bc in (True,False):
        B=S.build_blocks_z(N=N,zmin=zmin,zmax=zmax,edge_bc=bc); B0=S.build_blocks_z(N=N,zmin=zmin,zmax=zmax,edge_bc=bc,edge_tail=False)
        zz=B['z']; Fr=0.5*B['Q_op']; F0=0.5*B0['Q_op']
        Ff=B['M_GG']+2*B['Dz']       # forcing as used in M_GG
        Fx=0.5*(Qx[1:,1:] if bc else Qx)
        for g in gam:
            lr,vr=lam(Fr,zz,g); lf,_=lam(Ff,zz,g); lx,vx=lam(Fx,zz,g); l0,_=lam(F0,zz,g)
            out.append(f"N={N:4d} zmin={zmin:5.0f} zmax={zmax:3.0f} bc={int(bc)} g={g:.2f}  repo_tail(Q_op/2)={lr:+11.4e} repo_tail(forcing)={lf:+11.4e} EXACTKERNEL={lx:+11.4e} truncated={l0:+10.3f} | argmax|v| z: repo {zz[np.abs(vr).argmax()]:+.2f} exact {zz[np.abs(vx).argmax()]:+.2f}")
    return out
if __name__=="__main__":
    cfg=[(400,-10.,12.),(800,-10.,12.),(1600,-10.,12.),(1200,-12.,12.),(600,-12.,12.),(1600,-14.,12.),(800,-14.,12.),(800,-6.,12.),(800,-10.,7.),(1400,-16.,12.)]
    with Pool(4) as p:
        for o in p.imap(run,cfg):
            print("\n".join(o)); print(); sys.stdout.flush()
