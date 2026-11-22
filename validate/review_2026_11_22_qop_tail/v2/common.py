import sys, numpy as np
from scipy.integrate import quad
REPO="/Volumes/4TB SD/ClaudeCode/singular-blowup-1d"
sys.path.insert(0, REPO+"/src")
import stability_edge_z as S
import importlib.util
def load_old():
    sp=importlib.util.spec_from_file_location("old_sez","/private/tmp/qop_review_v2/old/src/stability_edge_z.py")
    m=importlib.util.module_from_spec(sp); sp.loader.exec_module(m); return m

def Ga(z): 
    w=np.exp(z); return w/(2+w)**2
def mk_bump(c,s=0.4):
    return lambda z: np.exp(-((z-c)/s)**2)
def mk_jump(zmin):
    return lambda z: np.exp(zmin-z)
def mk_azero(zmin):
    # two bumps with zero edge functional a = -(1/pi) int e^{-z/2} G
    b1=mk_bump(zmin+1.0,0.4); b2=mk_bump(zmin+3.0,0.4)
    i1=quad(lambda t: np.exp(-t/2)*b1(t), zmin-6, zmin+12, points=[zmin+1], limit=200)[0]
    i2=quad(lambda t: np.exp(-t/2)*b2(t), zmin-6, zmin+12, points=[zmin+3], limit=200)[0]
    return lambda z: b1(z)-(i1/i2)*b2(z)

def J_exact(G, z, lo, hi, extra_pts=()):
    """(1/pi) int_lo^hi e^{z'/2} G(z') ln|e^{z-z'}-1| dz'  = int_0^{e^z} H[g](1+y) dy, g supported [lo,hi]"""
    def f(zp):
        d=z-zp
        if d>0:
            L = d+np.log1p(-np.exp(-d)) if d>1e-300 else -700.
        elif d<0:
            L = np.log1p(-np.exp(d))
        else:
            return 0.0
        return np.exp(zp/2)*G(zp)*L
    tot=0.0
    brk=sorted(set([lo,hi]+[p for p in list(extra_pts)+[z] if lo<p<hi]))
    # refine around the singular point
    if lo<z<hi:
        for e in (1e-3,1e-1,1.0):
            for q in (z-e,z+e):
                if lo<q<hi: brk.append(q)
    elif z==lo:
        for e in (1e-3,1e-1,1.0):
            if lo+e<hi: brk.append(lo+e)
    brk=sorted(set(brk))
    for a,b in zip(brk[:-1],brk[1:]):
        if b-a>4:
            n=int(np.ceil((b-a)/4)); ed=np.linspace(a,b,n+1)
        else: ed=[a,b]
        for a2,b2 in zip(ed[:-1],ed[1:]):
            tot+=quad(f,a2,b2,limit=200,epsabs=1e-14,epsrel=1e-12)[0]
    return tot/np.pi
def Q_exact(G,z,lo,hi,extra_pts=()):
    return np.exp(-z)*J_exact(G,z,lo,hi,extra_pts)
