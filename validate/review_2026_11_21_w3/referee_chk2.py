import numpy as np, itertools
from scipy.special import loggamma
from scipy.linalg import eigh
t=np.linspace(-60,60,48001); dt=t[1]-t[0]
def logN(s): 
    # log of cos(pi s)Gamma(2s+2)/pi, stable
    lc=np.log(np.cos(np.pi*s.real)*np.cosh(np.pi*s.imag)-1j*np.sin(np.pi*s.real)*np.sinh(np.pi*s.imag)) if np.abs(s.imag).max()<200 else None
    return lc+loggamma(2*s+2)-np.log(np.pi)
def fam(sig,tau,ab):
    s=sig+1j*t; x=1-np.exp(-tau); lN=logN(s)
    cols=[np.exp(lN+loggamma(2*s+a)-(2*s+a)*np.log(b+x)) for a,b in ab]
    return np.array(cols).T
ab=[(a,b) for a in (0.5,1,2,3) for b in (0.25,0.5,1,2,4,8)]
for sig in (0.1,0.2,0.25,0.26,0.3,0.4,0.45):
    V0=fam(sig,0.0,ab); P0=(V0.conj().T@V0)*dt
    w,U=np.linalg.eigh(P0); keep=w>w.max()*1e-9; W=U[:,keep]/np.sqrt(w[keep])
    out=[]
    for tau in (0.02,0.1,0.3,1.0,3.0):
        V=fam(sig,tau,ab); P=(V.conj().T@V)*dt
        lam=np.linalg.eigvalsh(W.conj().T@P@W).max()
        out.append(np.sqrt(lam))
    print(f"sigma={sig}: rank {keep.sum()}  sup ||S G||e^(2 sig tau)/||G|| over span, tau=.02,.1,.3,1,3:", " ".join(f"{o:.5f}" for o in out))
