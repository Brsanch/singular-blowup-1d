import numpy as np
_x,_w=np.polynomial.legendre.leggauss(16)
def _kern(zi,zp):
    d=zi-zp
    with np.errstate(all='ignore'):
        e=-np.expm1(-np.abs(d))
        L=np.where(d>0,d,0.0)+np.log(e)
    L=np.where(np.isfinite(L),L,0.0)
    return np.exp(zp/2)*L
def exact_Q(z, chunk=100):
    """Collocation matrix of Q on hat functions, grid-supported g: 
       (Q G)(z_i)= e^{-z_i}/pi * int_{z0}^{zN} e^{z'/2} G(z') ln|e^{z_i-z'}-1| dz'  (independent of repo code)."""
    N=len(z); a=z[:-1]; b=z[1:]; h=b-a
    Q=np.zeros((N,N))
    zq=0.5*(a+b)[:,None]+0.5*h[:,None]*_x[None,:]      # (N-1,16)
    wq=0.5*h[:,None]*_w[None,:]
    lam=(zq-a[:,None])/h[:,None]
    for s in range(0,N,chunk):
        idx=np.arange(s,min(N,s+chunk))
        K=_kern(z[idx][:,None,None],zq[None,:,:])*wq[None]       # (m,N-1,16)
        WL=(K*(1-lam)[None]).sum(-1); WR=(K*lam[None]).sum(-1)   # (m,N-1)
        # zero out adjacent cells (singular), redo below
        for r,i in enumerate(idx):
            for c in (i-1,i):
                if 0<=c<N-1:
                    WL[r,c]=0;WR[r,c]=0
                    # graded panels toward z_i inside cell c
                    if c==i:   # singular at left end a_c
                        ed=a[c]+h[c]*np.concatenate(([0.0],0.3**np.arange(34,-1,-1)))
                    else:      # singular at right end
                        ed=b[c]-h[c]*np.concatenate(([0.0],0.3**np.arange(34,-1,-1)))[::-1]
                    pa=ed[:-1];pb=ed[1:]
                    zz=0.5*(pa+pb)[:,None]+0.5*(pb-pa)[:,None]*_x[None,:]
                    ww=0.5*(pb-pa)[:,None]*_w[None,:]
                    kk=_kern(z[i],zz)*ww
                    ll=(zz-a[c])/h[c]
                    WL[r,c]=(kk*(1-ll)).sum(); WR[r,c]=(kk*ll).sum()
        Q[idx,:-1]+=WL; Q[idx,1:]+=WR
    return np.exp(-z)[:,None]*Q/np.pi
if __name__=="__main__":
    from common import *; from fastq import Q_fast
    zmin=-10.;z=np.linspace(zmin,7.,400)
    Qx=exact_Q(z)
    for nm,G in (('Ga',Ga),('bump',mk_bump(zmin+1)),('jump',mk_jump(zmin))):
        ref=Q_fast(G,z,zmin,7.); q=Qx@G(z)
        print(nm, np.abs(q-ref).max()/np.abs(ref).max())
    z=np.linspace(zmin,7.,1600); Qx=exact_Q(z)
    for nm,G in (('Ga',Ga),('bump',mk_bump(zmin+1)),('jump',mk_jump(zmin))):
        ref=Q_fast(G,z,zmin,7.); q=Qx@G(z)
        print(nm, np.abs(q-ref).max()/np.abs(ref).max())
