import numpy as np
from scipy.special import loggamma
def N(s): return np.cos(np.pi*s)*np.exp(loggamma(2*s+2))/np.pi
def c(w): return -np.tan(np.pi*w)/(2*w+1)
t=np.linspace(-12,12,4801); dt=t[1]-t[0]
z=np.linspace(-60,60,60001); dz=z[1]-z[0]
def run(sig,s0,L):
    # form identity at tau=0
    Gh=lambda s: N(s)*np.exp(L*(s-s0)**2)
    s=sig+1j*t
    form_a=np.real(np.sum(np.conj(Gh(s))*c(s+0.5)*Gh(s+0.5)))*dt
    w=sig+0.25+1j*t
    form_b=np.sum(np.real(c(w))*np.abs(Gh(w))**2)*dt
    nrm=np.sum(np.abs(Gh(s))**2)*dt
    # real-space evolution
    k0=lambda zz: np.exp(s0*zz-zz**2/(4*L))/np.sqrt(4*np.pi*L)
    E=np.exp(-np.outer(sig+1j*t,z))  # big: 4801 x 60001 complex -> too big; chunk
    return form_a,form_b,nrm
def evol(sig,s0,L,taus):
    k0=lambda zz: np.exp(s0*zz-zz**2/(4*L))/np.sqrt(4*np.pi*L)
    s=sig+1j*t; Ns=N(s); out=[]
    for tau in taus:
        x=1-np.exp(-tau); zz=z-2*tau
        k=k0(zz)*np.exp(-x*np.exp(-zz/2))
        g=np.exp(-sig*z)*k
        kh=np.array([np.sum(np.exp(-1j*tt*z)*g)*dz for tt in t[::4]])
        out.append(np.sqrt(np.sum(np.abs(Ns[::4]*kh)**2))*np.exp(2*sig*tau))
    return np.array(out)
taus=[0,0.01,0.03,0.1,0.3,1.0]
for sig,s0,L in [(0.4,0.65,4.0),(0.4,0.65,20.0),(0.3,0.55,20.0),(0.45,0.7,20.0),(0.2,0.45,20.0),(0.2,0.3,4.0)]:
    fa,fb,n=run(sig,s0,L)
    e=evol(sig,s0,L,taus)
    print(f"sig={sig} s0={s0} L={L}: form(line sig)={fa/n:+.5f} form(line sig+1/4)={fb/n:+.5f} | ||S G||e^(2 sig tau)/||G|| at tau={taus[1:]}: "+" ".join(f"{v:.5f}" for v in e[1:]/e[0]),
          f"| check tau=0 norm ratio {e[0]**2*4*dt/n:.6f}")
