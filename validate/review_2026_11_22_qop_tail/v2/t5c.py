import numpy as np, sys
from multiprocessing import Pool
from common import *; from exactmat import exact_Q; from fastq import Q_fast
def lam(F,z,g):
    w=np.exp(-g*z); A=(w[:,None]*F)/w[None,:]; return np.linalg.eigvalsh(0.5*(A+A.T))[-1]
def ramp_A(Q,z,dz,gamma=0.4,z0=-6.,s=2.):
    pp=2*gamma*0.5*(1+np.tanh((z-z0)/s)); phi=np.cumsum(pp)*dz; sw=np.exp(-0.5*phi)
    Ft=(sw[:,None]*(0.5*Q))/sw[None,:]
    return np.linalg.eigvalsh(-np.diag(pp)+0.5*(Ft+Ft.T))[-1]
def corrected_Q(N,zmin,zmax):
    """repo construction but with the Hilbert diagonal corrected to a grid-supported g
       (remove the implicit 'g = g_i on [1,Y_0]' term), free-edge matrices, then drop node 0."""
    z=np.linspace(zmin,zmax,N); dz=z[1]-z[0]; Y=1+np.exp(z); ym1=Y-1; Ymax=Y[-1]+0.5*(Y[-1]-Y[-2])
    Hz=S.hilbert_support_z(Y,Ymax)
    corr=np.zeros(N); corr[1:]=np.log(ym1[1:]/(ym1[1:]-ym1[0]))/np.pi
    Hz=Hz-np.diag(corr)
    Hg=Hz@np.diag(ym1**-0.5)
    EA=np.zeros((N,N))
    for i in range(1,N):
        tw=np.full(i+1,dz);tw[0]*=.5;tw[i]*=.5; EA[i,:i+1]=np.exp(z[:i+1]-z[i])*tw
    Q=EA@Hg+np.outer(np.exp(-z),S.edge_tail_row(z))
    return Q[1:,1:],z[1:]
def run(a):
    N,zmin=a; zmax=12.
    B=S.build_blocks_z(N=N,zmin=zmin,zmax=zmax,edge_bc=True); z=B['z']
    Qc,_=corrected_Q(N,zmin,zmax)
    s=f"N={N} zmin={zmin}: ramp A_full repo={ramp_A(B['Q_op'],z,B['dz']):+.4f} corrected-H={ramp_A(Qc,z,B['dz']):+.4f}"
    if N<=2000:
        Qx=exact_Q(np.linspace(zmin,zmax,N))[1:,1:]
        s+=f" exactkernel={ramp_A(Qx,z,B['dz']):+.4f}\n"
        for g in (0.25,0.30,0.40):
            s+=f"     g={g}: C_force repo={lam(.5*B['Q_op'],z,g):+.4e} corrected-H={lam(.5*Qc,z,g):+.4e} exact={lam(.5*Qx,z,g):+.4e}\n"
        G=mk_bump(zmin+1); ref=Q_fast(G,z,zmin,zmax)
        s+=f"     bump@zmin+1 rel Q err: repo={np.abs(B['Q_op']@G(z)-ref).max()/np.abs(ref).max():.2e} corrected-H={np.abs(Qc@G(z)-ref).max()/np.abs(ref).max():.2e} exact={np.abs(Qx@G(z)-ref).max()/np.abs(ref).max():.2e}"
    return s
if __name__=="__main__":
    cfg=[(700,-10.),(1400,-10.),(2000,-10.),(2800,-10.),(4000,-10.),(1400,-8.),(2800,-8.),(4000,-12.),(1400,-12.)]
    with Pool(4) as p:
        for o in p.imap(run,cfg): print(o); sys.stdout.flush()
