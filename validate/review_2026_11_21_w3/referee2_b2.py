import numpy as np
from scipy.special import loggamma
from scipy.linalg import eigh
np.seterr(all='ignore')
def N(s): return np.cos(np.pi*s)*np.exp(loggamma(2*s+2))/np.pi
def gen_max(A,B,tol):
    w,V=eigh(B); keep=w>tol*w.max(); P=V[:,keep]/np.sqrt(w[keep])
    lam,U=eigh(P.conj().T@A@P); return lam[-1], keep.sum(), P@U[:,-1]
def run(sig,params,taus,T=250,n=50001,tol=1e-8):
    t=np.linspace(-T,T,n); dt=t[1]-t[0]
    s=sig+1j*t; sg=np.where(t>=0,1.0,-1.0); lN=-1j*np.pi*s*sg+np.log((1+np.exp(2j*np.pi*s*sg))/2)+loggamma(2*s+2)-np.log(np.pi)
    def lk(a,b,tau,sh=0): 
        ss=s+sh; return -2*ss*tau+np.log(2)+loggamma(2*ss+a)-(2*ss+a)*np.log(b+1-np.exp(-tau))
    G0=np.array([np.exp(lN+lk(a,b,0)) for a,b in params])
    nr=np.sqrt(np.sum(abs(G0)**2,axis=1)*dt/(2*np.pi)); G0/=nr[:,None]
    tail=np.max(abs(G0[:,[0,-1]]))
    B=(G0.conj()@G0.T)*dt/(2*np.pi); B=(B+B.conj().T)/2
    MG=np.array([np.exp(lN+lk(a,b,0))*(-2*s+2*sig)-np.exp(lN+lk(a,b,0,0.5)) for a,b in params])/nr[:,None]
    F=(G0.conj()@MG.T)*dt/(2*np.pi); F=(F+F.conj().T)/2
    fm,r,_=gen_max(F,B,tol)
    best=(0,None)
    for tau in taus:
        Gt=np.array([np.exp(lN+lk(a,b,tau)) for a,b in params])/nr[:,None]
        A=(Gt.conj()@Gt.T)*dt/(2*np.pi); A=(A+A.conj().T)/2
        lam,_,v=gen_max(A,B,tol)
        # recompute ratio directly from vector (avoid Gram roundoff)
        g0=v@G0; gt=v@Gt
        rat=np.sqrt(np.sum(abs(gt)**2)/np.sum(abs(g0)**2))*np.exp(2*sig*tau)
        if rat>best[0]: best=(rat,tau)
    return fm,r,best,tail
fams={
 'real':[(a,b) for a in [0.2,0.5,1,2,3,4.5,6] for b in [0.05,0.3,1,4,20]],
 'complex b':[(a,r*np.exp(1j*th)) for a in [0.3,1,2.5,5] for r in [0.3,1,3] for th in [-1.2,-0.6,0,0.6,1.2]],
 'small a':[(a,b) for a in [0.02,0.05,0.1,0.2,0.4] for b in [0.1,1,10,100]],
 'large a':[(a,b) for a in [6,8,10,12] for b in [0.5,1,2,4]],
 'near-imag b':[(a,np.exp(1j*th)) for a in [0.1,0.5,1,2] for th in [-1.45,-1.3,1.3,1.45]],
}
taus=[1e-3,1e-2,0.05,0.2,0.5,1,2,4]
for sig in [0.02,0.1,0.2,0.25,0.26,0.27,0.3,0.4]:
    for nm,p in fams.items():
        T=250 if nm!='near-imag b' else 600
        fm,r,best,tail=run(sig,p,taus,T=T,n=int(T*200)+1)
        print(f"sig={sig:4.2f} {nm:12s} rank={r:2d} maxform={fm:+.3e} max ratio={best[0]:.9f} at tau={best[1]} tail={tail:.1e}")
