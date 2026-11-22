import numpy as np, sys
from common import *
from scipy.integrate import quad, dblquad
# (i) kernel formula vs truly independent double integral int_0^{vmin} H[g](1+y) dy
def tail_double(G, zmin, zmax, pts=()):
    vmin=np.exp(zmin)
    def Hg(y):  # (1/pi) int g(v')/(y-v') dv', v'=e^{z'}, g=e^{-z'/2}G
        f=lambda zp: np.exp(zp/2)*G(zp)/(y-np.exp(zp))
        brk=sorted(set([zmin,zmax]+[p for p in pts if zmin<p<zmax]+[zmin+e for e in (1e-6,1e-4,1e-2,0.3,2.,6.) if zmin+e<zmax]))
        return sum(quad(f,a,b,limit=400,epsabs=0,epsrel=1e-11)[0] for a,b in zip(brk[:-1],brk[1:]))/np.pi
    # y in [0,vmin], log-singular at y->vmin if G(zmin)!=0: substitute y = vmin(1-e^{-s})
    f2=lambda s: Hg(vmin*(1-np.exp(-s)))*vmin*np.exp(-s)
    return sum(quad(f2,a,b,limit=200,epsrel=1e-10)[0] for a,b in [(0,1),(1,5),(5,12),(12,30)])
def tail_kernel(G,zmin,zmax,pts=()):
    return J_exact(G,zmin,zmin,zmax,pts)
for zmin in (-10.,-14.):
    zmax=7.
    fs={'a w/(2+w)^2':(Ga,()), 'b bump@zmin+1':(mk_bump(zmin+1),(zmin+1,)), 'c jump e^{zmin-z}':(mk_jump(zmin),()), 'd a=0 pair':(mk_azero(zmin),(zmin+1,zmin+3))}
    for name,(G,pts) in fs.items():
        td=tail_double(G,zmin,zmax,pts); tk=tail_kernel(G,zmin,zmax,pts)
        a_fun=-(1/np.pi)*quad(lambda t:np.exp(-t/2)*G(t),zmin,zmax,points=[p for p in pts],limit=400)[0]
        print(f"zmin={zmin} {name:20s} double={td:+.10e} kernelquad={tk:+.10e} rel={abs(td-tk)/abs(td):.1e}  hand a*vmin={a_fun*np.exp(zmin):+.4e} ratio hand/exact={a_fun*np.exp(zmin)/tk:.4f}")
        for N in (200,400,800,1600,3200):
            z=np.linspace(zmin,zmax,N); t=S.edge_tail_row(z)
            v=t@G(z); v1=t[1:]@G(z[1:])
            print(f"     N={N:5d} row@G={v:+.10e} relerr={abs(v-tk)/abs(tk):.2e}   t[0]={t[0]:+.3e} t[0]G0/total={t[0]*G(z[0])/v:+.3e}")
