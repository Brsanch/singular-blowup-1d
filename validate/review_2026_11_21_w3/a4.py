# Attack 4: sup over a large explicit family of ||S(tau)g||_sigma/||g||_sigma.
# g_j^ (s) = N(s) 2 Gamma(2s+a) b^{-2s-a};  S(tau): b -> b+x, times e^{-2 s tau}, x=1-e^{-tau}.
# (a) first verify pointwise that this closed form solves d/dtau G^(s) = -2s G^(s) + c(s+1/2) G^(s+1/2)  (no conjugation used)
import numpy as np, sys
from scipy.special import loggamma
from scipy.linalg import qr, svd
def logN(s):
    t=s.imag; pos=t>=0
    e=np.where(pos,np.exp(2j*np.pi*np.where(pos,s,1j)),np.exp(-2j*np.pi*np.where(pos,-1j,s)))
    lc=np.where(pos,-1j*np.pi*s,1j*np.pi*s)-np.log(2)+np.log1p(e)
    return lc+loggamma(2*s+2)-np.log(np.pi)
def logG(s,a,b,tau):
    x=1-np.exp(-tau)
    return logN(s)+np.log(2)+loggamma(2*s+a)-2*s*tau-(2*s+a)*np.log(b+x)
c=lambda s:-np.tan(np.pi*s)/(2*s+1)
# PDE check
worst=0
for a,b in [(1.0,1.0),(0.5,0.3+0.8j),(2.0,0.05-0.2j),(-0.3,3+1j)]:
    for s in np.array([0.3+0.7j,0.3-3j,0.3+0.01j,0.45+2j]):
        for tau in [0.01,0.5,2.0]:
            h=1e-5
            d=(np.exp(logG(np.array([s]),a,b,tau+h))-np.exp(logG(np.array([s]),a,b,tau-h)))/(2*h)
            rhs=-2*s*np.exp(logG(np.array([s]),a,b,tau))+c(s+0.5)*np.exp(logG(np.array([s+0.5]),a,b,tau))
            worst=max(worst,abs(d-rhs)[0]/abs(rhs)[0])
print("transform-side evolution eq. residual (FD in tau), worst rel:",worst)

def run(sig,avals,rs,betas,taus,T=500.,dt=0.02,thr=1e-11,verbose=False):
    t=np.arange(-T,T+dt/2,dt); s=sig+1j*t
    params=[(a,r*np.exp(1j*be)) for a in avals for r in rs for be in betas]
    def mat(tau):
        cols=[]
        for a,b in params:
            cols.append(np.exp(logG(s,a,b,tau)))
        return np.array(cols).T*np.sqrt(dt)
    P0=mat(0.0)
    nrm=np.linalg.norm(P0,axis=0)
    P0n=P0/nrm
    U,S,Vh=svd(P0n,full_matrices=False)
    keep=S>thr*S[0]
    B=Vh[keep].conj().T/S[keep]      # coefficients c = B y, ||P0n c|| = ||y||
    out=[]
    for tau in taus:
        Pt=mat(tau)/nrm
        A=Pt@B
        sv=svd(A,compute_uv=False)
        single=np.linalg.norm(Pt,axis=0).max()
        out.append((tau,sv[0],sv[0]*np.exp(2*sig*tau),single))
    return len(params),int(keep.sum()),out
if __name__=="__main__" and False:
    taus=[1e-4,1e-3,1e-2,0.03,0.1,0.2,0.4,0.7,1.0,1.5,2.5,4.0,8.0]
    for sig in [0.05,0.2,0.3,0.4,0.45,0.49]:
        for label,avals,rs,betas in [
            ("small",[1.0],np.exp(np.linspace(-2,2,5)),[0.0,0.8,-0.8]),
            ("medium",[1.0],np.exp(np.linspace(-6,4,11)),[0,0.6,-0.6,1.2,-1.2,1.45,-1.45]),
            ("large",[1.0,0.3,2.5],np.exp(np.linspace(-10,5,16)),[0,0.6,-0.6,1.2,-1.2,1.45,-1.45,1.53,-1.53]),
        ]:
            n,k,out=run(sig,avals,rs,betas,taus)
            best=max(out,key=lambda o:o[1]); bestw=max(out,key=lambda o:o[2])
            print(f"sigma={sig} family={label} n={n} rank_kept={k}: max ratio={best[1]:.4f} at tau={best[0]}; max ratio*e^(2 sig tau)={bestw[2]:.4f} at tau={bestw[0]}")
            if label=="large":
                print("   tau: ratio, ratio*e^{2 sig tau}: "+"; ".join(f"{o[0]:g}: {o[1]:.4f}, {o[2]:.4f}" for o in out))
            sys.stdout.flush()
