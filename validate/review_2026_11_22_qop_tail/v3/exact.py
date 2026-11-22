import numpy as np
from scipy.integrate import quad
def phi(t):
    t=np.asarray(t,float); out=np.ones_like(t); m=abs(t)>1e-8
    out[m]=t[m]/np.expm1(t[m]); return out
def exact_Hg(G,z,a,b,brk=()):
    """(1/pi) pv int_a^b e^{z'/2}G(z')/(e^z-e^{z'}) dz'"""
    h=lambda zp: float(-np.exp(0.5*zp-z)*G(zp)*phi(np.array([zp-z]))[0])
    pts=sorted(set([a,b]+[p for p in brk if a<p<b]))
    # piece containing z uses cauchy weight
    tot=0.0
    for lo,hi in zip(pts[:-1],pts[1:]):
        if lo<z<hi:
            d=min(z-lo,hi-z,0.5)
            tot+=quad(h,z-d,z+d,weight='cauchy',wvar=z,limit=400,epsabs=1e-13,epsrel=1e-12)[0]
            f=lambda zp:h(zp)/(zp-z)
            if z-d>lo: tot+=quad(f,lo,z-d,limit=400,epsabs=1e-13,epsrel=1e-12)[0]
            if z+d<hi: tot+=quad(f,z+d,hi,limit=400,epsabs=1e-13,epsrel=1e-12)[0]
        else:
            f=lambda zp:h(zp)/(zp-z)
            tot+=quad(f,lo,hi,limit=400,epsabs=1e-13,epsrel=1e-12)[0]
    return tot/np.pi
def exact_Q(G,z,a,b,brk=()):
    f=lambda zp: np.exp(0.5*zp)*G(zp)*np.log(abs(np.expm1(z-zp)))
    pts=sorted(set([a,b]+[p for p in list(brk)+[z] if a<p<b]))
    tot=0.0
    for lo,hi in zip(pts[:-1],pts[1:]):
        tot+=quad(f,lo,hi,limit=400,epsabs=1e-14,epsrel=1e-12)[0]
    return np.exp(-z)*tot/np.pi
def bump(c,r):
    def B(zp):
        x=(np.asarray(zp,float)-c)/r
        out=np.zeros_like(x); m=abs(x)<1
        out[m]=np.exp(1-1/(1-x[m]**2)); return out if out.ndim else float(out)
    return B
