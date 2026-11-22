import numpy as np
from scipy.special import loggamma
from scipy.linalg import eigh
np.seterr(all='ignore')
def N(s): return np.cos(np.pi*s)*np.exp(loggamma(2*s+2))/np.pi
def c(w): return -np.tan(np.pi*w)/(2*w+1)
def rec(x,t): return -((2*x+1)*np.sin(2*np.pi*x)+2*t*np.sinh(2*np.pi*t))/((np.cos(2*np.pi*x)+np.cosh(2*np.pi*t))*abs(2*(x+1j*t)+1)**2)
# formula check
rng=np.random.default_rng(0)
xs=rng.uniform(-0.4,0.95,2000); ts=rng.uniform(-3,3,2000)
print("Re c formula max abs err", np.max(np.abs(c(xs+1j*ts).real-rec(xs,ts))/(1+np.abs(rec(xs,ts)))))
# sign checks
X,T=np.meshgrid(np.linspace(1e-4,0.5-1e-4,400),np.linspace(-6,6,1201))
print("max Re c on 0<x<1/2:", rec(X,T).max())
t=np.linspace(-6,6,1200); print("max Re c on x=1/2 (t!=0):", rec(0.5,t).max(), " limit t->0:", rec(0.5,1e-5), -1/(2*np.pi))
for x in [0.51,0.6,0.7,0.749]: print("x",x,"Re c(x)=",rec(x,0.0), "Re c(x+2i)=",rec(x,2.0))
# data families: khat(s)
def gam(a,b): return lambda s: 2*np.exp(loggamma(2*s+a)-(2*s+a)*np.log(b))
def gauss(z0,w,om): return lambda s: np.sqrt(2*np.pi)*w*np.exp(-(s-1j*om)*z0+w*w*(s-1j*om)**2/2)
def combo(fs,cs): return lambda s: sum(ci*f(s) for ci,f in zip(cs,fs))
t=np.linspace(-60,60,240001); dt=t[1]-t[0]
def two_lines(kh,sig):
    s=sig+1j*t
    G=N(s)*kh(s); MG=N(s)*(-2*s*kh(s)-kh(s+0.5))
    lhs=np.sum(((MG+2*sig*G)*np.conj(G)).real)*dt/(2*np.pi)
    nrm=np.sum(abs(G)**2)*dt/(2*np.pi)
    w=sig+0.25+1j*t
    Gw=N(w)*kh(w); rc=rec(sig+0.25,t)
    integ=rc*abs(Gw)**2
    integ[~np.isfinite(integ)]=0
    rhs=np.sum(integ)*dt/(2*np.pi)
    return lhs,rhs,nrm
data={'gam(1,1)':gam(1,1),'gam(0.3,2+1.5j)':gam(0.3,2+1.5j),'gam(4,0.2-0.1j)':gam(4,0.2-0.1j),
 'gauss(0,1,0)':gauss(0,1,0),'gauss(-3,0.5,2)':gauss(-3,0.5,2),'gauss(2,2,-1)':gauss(2,2,-1),'gauss(0,3,0.0)':gauss(0,3,0),
 'mix':combo([gam(1,1),gam(2.5,3-2j),gauss(-1,0.7,1.3)],[1,-0.7+0.2j,0.4j])}
print("\n(a) identity two lines: sigma, data, LHS/||G||^2, RHS/||G||^2, relerr")
for sig in [0.03,0.1,0.2,0.25,0.27,0.3,0.4,0.47]:
    for nm,kh in data.items():
        l,r,n=two_lines(kh,sig); print(f"{sig:5.2f} {nm:18s} {l/n:+.10e} {r/n:+.10e} {abs(l-r)/abs(n):.1e}")
