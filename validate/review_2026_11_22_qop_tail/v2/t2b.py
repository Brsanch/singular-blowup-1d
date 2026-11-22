import numpy as np
from common import *; from fastq import *
from scipy.integrate import quad
def H_exact(G,z,lo,hi,pts=()):
    # (1/pi) pv int_lo^hi e^{z'/2}G(z')/(e^z - e^{z'}) dz'  for lo<z<hi
    v=np.exp(z); gz=np.exp(z/2)*G(z)
    f=lambda zp: (np.exp(zp/2)*G(zp)-gz*np.exp(zp-z))/(v-np.exp(zp)) if zp!=z else 0.0   # subtract g(z)-like term: e^{z'/2}G -> use h(z')=e^{zp/2}G(zp); pv of h(z)e^{zp-z}/(v-e^{zp}) analytic
    brk=sorted(set([lo,hi,z]+[p for p in pts if lo<p<hi]+[q for e in (1e-2,0.3,2.) for q in (z-e,z+e) if lo<q<hi]))
    I=sum(quad(f,a,b,limit=300,epsrel=1e-11)[0] for a,b in zip(brk[:-1],brk[1:]))
    # pv int_lo^hi e^{zp-z}/(v-e^{zp}) dz' = (1/v) pv int_{vlo}^{vhi} dv'/(v-v') = (1/v) ln((v-vlo)/(vhi-v))
    I+=gz/v*np.log((v-np.exp(lo))/(np.exp(hi)-v))
    return I/np.pi
for zmin in (-10.,):
  for name,G,pts in (('bump',mk_bump(zmin+1),(zmin+1,)),('Ga',Ga,())):
    for N in (800,3200):
        B=S.build_blocks_z(N=N,zmin=zmin,zmax=7.,edge_bc=False); z=B['z']
        hg=B['Hg_op']@G(z)
        idx=[int(round((d)/(z[1]-z[0]))) for d in (0.25,0.5,1.0,1.2,2.0,3.0,5.0)]
        print(name,"N",N)
        for i in idx:
            he=H_exact(G,z[i],zmin,7.,pts); hc=H_exact(G,z[i],-60.,7.,pts)
            gi=np.exp(-z[i]/2)*G(z[i]); pred=gi*np.log(np.exp(z[i])/(np.exp(z[i])-np.exp(zmin)))/np.pi
            print(f"  z-zmin={z[i]-zmin:5.2f} Hg_op={hg[i]:+.5e} exact(grid-supp)={he:+.5e} exact(cont)={hc:+.5e} diff={hg[i]-he:+.4e} predicted g_i ln(v_i/(v_i-v0))/pi={pred:+.4e}")
