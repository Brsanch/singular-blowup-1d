import numpy as np
from scipy.special import loggamma
np.seterr(all='ignore')
def absN2(x,t):
    s=x+1j*t; sg=np.where(t>=0,1.0,-1.0)
    lN=-1j*np.pi*s*sg+np.log((1+np.exp(2j*np.pi*s*sg))/2)+loggamma(2*s+2)-np.log(np.pi)
    return np.exp(2*lN.real)
def rec(x,t):
    # stable for large t
    e=np.exp(-2*np.pi*np.abs(t))
    num=(2*x+1)*np.sin(2*np.pi*x)*2*e+2*t*np.sign(t)*(1-e*e)
    den=(np.cos(2*np.pi*x)*2*e+1+e*e)*abs(2*(x+1j*t)+1)**2
    return -num/den
def bump(y): 
    out=np.zeros_like(y); m=abs(y)<1; out[m]=np.exp(-1/(1-y[m]**2)); return out
def FT(f,y):  # returns t, int f e^{-ity} dy
    n=len(y); dy=y[1]-y[0]
    t=2*np.pi*np.fft.fftfreq(n,dy)
    F=np.fft.fft(f)*dy*np.exp(-1j*t*y[0])
    return t,F
def test(sig,R,om=0.0,n=2**17,taus=(1e-5,1e-4,1e-3,1e-2,0.1,0.3,1,3)):
    x=sig+0.25; L=4*R
    y=np.linspace(-L,L,n,endpoint=False)
    eta=bump(y/R)*np.exp(1j*om*y)          # = e^{-x y} k_R
    k_sig=eta*np.exp(0.25*y)               # = e^{-sig y} k_R
    t,K=FT(k_sig,y); dtt=abs(t[1]-t[0])
    w=absN2(sig,t); nrm=np.sum(w*abs(K)**2)*dtt/(2*np.pi)
    # form via line sigma+1/4 (item 1 RHS)
    t2,K2=FT(eta,y); form=np.sum(rec(x,t2)*absN2(x,t2)*abs(K2)**2)*dtt/(2*np.pi)
    # form via line sigma (LHS):  -Re int N(s) khat(s+1/2) conj(N khat)
    # khat(s+1/2) on Re s = sig : FT of e^{-(sig+1/2)y}k = eta*e^{-y/4}
    _,K3=FT(eta*np.exp(-0.25*y),y)
    s=sig+1j*t; sg=np.where(t>=0,1.0,-1.0)
    lN=-1j*np.pi*s*sg+np.log((1+np.exp(2j*np.pi*s*sg))/2)+loggamma(2*s+2)-np.log(np.pi)
    Ns=np.exp(lN)
    lhs=np.sum((-(Ns*K3)*np.conj(Ns*K)).real)*dtt/(2*np.pi)
    res=[]
    for tau in taus:
        xx=1-np.exp(-tau)
        _,Kt=FT(k_sig*np.exp(-xx*np.exp(-y/2)),y)
        r2=np.sum(w*abs(Kt)**2)/np.sum(w*abs(K)**2)   # = e^{4 sig tau}||S G||^2/||G||^2
        res.append((tau,np.sqrt(r2),(r2-1)/tau))
    print(f"sig={sig} R={R} om={om}: form/||G||^2 line(sig+1/4)={form/nrm:+.6e} line(sig)={lhs/nrm:+.6e}; 2*form/||G||^2={2*form/nrm:+.6e}")
    print("    tau: ratio, (ratio^2-1)/tau :", "  ".join(f"{a:g}: {b:.9f}, {c_:+.4e}" for a,b,c_ in res))
for sig,R,om in [(0.27,10,0),(0.27,30,0),(0.27,60,0),(0.27,120,0),(0.26,120,0),(0.3,40,0),(0.4,20,0),
                 (0.25,60,0),(0.25,200,0),(0.25,20,0.5),(0.2,60,0),(0.1,30,0),(0.1,30,2.0),(0.03,60,0),(0.27,60,1.0)]:
    test(sig,R,om)
