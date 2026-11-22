"""Independent collocation of Q from (QG)(z)=(e^{-z}/pi) int e^{z'/2} G(z') ln|e^{z-z'}-1| dz',
G piecewise linear on the grid, zero outside [z_0, z_{N-1}].  No repo code used."""
import numpy as np
def lnabs_expm1(t):
    t=np.asarray(t,float); out=np.empty_like(t)
    p=t>0; out[p]=t[p]+np.log(-np.expm1(-t[p])); out[~p]=np.log(-np.expm1(t[~p])); return out
def smooth_part(t):   # ln|(e^t-1)/t|
    t=np.asarray(t,float); out=np.zeros_like(t); m=abs(t)>1e-6
    out[m]=lnabs_expm1(t[m])-np.log(abs(t[m])); out[~m]=t[~m]/2; return out
def Qmat(z,ng=10,chunk=200):
    N=len(z); xg,wg=np.polynomial.legendre.leggauss(ng)
    a=z[:-1]; b=z[1:]; h=b-a
    zq=a[:,None]+0.5*h[:,None]*(xg[None,:]+1)        # (N-1,ng)
    lam=(zq-a[:,None])/h[:,None]; wq=0.5*h[:,None]*wg[None,:]
    eq=np.exp(0.5*zq)
    Q=np.zeros((N,N))
    for s in range(0,N,chunk):
        rows=np.arange(s,min(N,s+chunk))
        t=z[rows][:,None,None]-zq[None,:,:]
        with np.errstate(divide='ignore'):
            k=lnabs_expm1(t)*eq[None]*wq[None]
        I0=(k*(1-lam)[None]).sum(-1); I1=(k*lam[None]).sum(-1)   # (rows, cells)
        # zero out the adjacent (singular) cells; redo them below
        for r_,i in enumerate(rows):
            for c in (i-1,i):
                if 0<=c<N-1: I0[r_,c]=0; I1[r_,c]=0
        Q[rows,:-1]+=I0; Q[rows,1:]+=I1
    # singular cells: ln|e^t-1| = ln|t| + smooth(t)
    xs,ws=np.polynomial.legendre.leggauss(40); u=0.5*(xs+1); wu=0.5*ws; p=6
    for i in range(N):
        for c,side in ((i-1,-1),(i,+1)):
            if not (0<=c<N-1): continue
            hc=h[c]
            # smooth part with GL
            zz=zq[c]; t=z[i]-zz; ks=smooth_part(t)*eq[c]*wq[c]
            Q[i,c]+=(ks*(1-lam[c])).sum(); Q[i,c+1]+=(ks*lam[c]).sum()
            # log part: s=|z'-z_i| in (0,hc), s=hc*u^p
            sdist=hc*u**p; ds=hc*p*u**(p-1)*wu
            zp=z[i]+side*sdist; lg=np.log(sdist)
            l=(zp-z[c])/hc
            kk=lg*np.exp(0.5*zp)*ds
            Q[i,c]+=(kk*(1-l)).sum(); Q[i,c+1]+=(kk*l).sum()
    return Q*np.exp(-z)[:,None]/np.pi
def grid(N,zmin,zmax): return np.linspace(zmin,zmax,N)
def csym(F,z,gam):
    sw=np.exp(-gam*z); At=(sw[:,None]*F)/sw[None,:]
    return float(np.linalg.eigvalsh(0.5*(At+At.T))[-1])
if __name__=="__main__":
    # self-test vs scipy quad on test functions
    from exact import exact_Q, bump
    for N in (200,400,800):
        zmin=-8.;z=grid(N,zmin,7.);Q=Qmat(z)
        for nm,G,a,b in (('bump',bump(zmin+1,.8),zmin+.2,zmin+1.8),('jump',lambda t:np.exp(zmin-np.asarray(t)),zmin,7.)):
            Gv=np.array([float(G(t)) for t in z]); idx=[0,1,3,N//20,N//10,N//4,N//2,N-2]
            # exact for the piecewise-linear interpolant: tests the quadrature itself
            Gi=lambda t:float(np.interp(t,z,Gv))
            ex=np.array([exact_Q(Gi,z[i],a,b,brk=list(z[(z>a)&(z<b)][::1]) if N<=200 else ()) for i in idx]) if N<=200 else None
            ex2=np.array([exact_Q(G,z[i],a,b) for i in idx])
            d=(Q@Gv)[idx]
            print(N,nm,"vs true G: max rel err",abs(d-ex2).max()/abs(ex2).max(), "" if ex is None else f"| vs interpolant (quadrature only): {abs(d-ex).max()/abs(ex).max():.2e}")
