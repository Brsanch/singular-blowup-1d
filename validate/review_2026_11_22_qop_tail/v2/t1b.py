import numpy as np
from common import *
from scipy.integrate import quad
import mpmath as mp
for zmin,zmax,N in ((-10.,7.,50),(-14.,12.,800),(-4.,3.,17)):
    z=np.linspace(zmin,zmax,N); t=S.edge_tail_row(z)
    k=lambda zp: (mp.e**(zp/2)*(-mp.log(-mp.expm1(zmin-zp)))) if zp>zmin else mp.mpf(0)
    mp.mp.dps=30
    pts=[zmin]+[zmin+10.0**(-j) for j in range(8,0,-1)]+list(np.arange(zmin+1,zmax,1.0))+[zmax]
    I0=float(-mp.quad(k,pts)/mp.pi); I1=float(-mp.quad(lambda x:k(x)*x,pts)/mp.pi)
    # single-weight checks: first, second, last weights vs mp
    h=z[1]-z[0]
    w0=float(-mp.quad(lambda x:k(x)*(z[1]-x)/h,[zmin,zmin+h*1e-6,zmin+h*1e-3,zmin+h*0.1,z[1]])/mp.pi)
    wl=float(-mp.quad(lambda x:k(x)*(x-z[-2])/h,[z[-2],z[-1]])/mp.pi)
    print(f"zmin={zmin} N={N}: sum(t) rel err={abs(t.sum()-I0)/abs(I0):.1e}  t@z rel err={abs(t@z-I1)/abs(I1):.1e}  t[0] rel err={abs(t[0]-w0)/abs(w0):.1e}  t[-1] rel err={abs(t[-1]-wl)/abs(wl):.1e}")
