import numpy as np, sys
from lib import *
from solver import *
sc = G0_maxnorm_scale()
def diag(zmin,zmax,n,nst,A,filt=None,S=1.5,verbose=False):
    s=Solver(zmin,zmax,n,filt=filt); h=A*h0_family(s)/sc; s.setup(S/nst)
    xi=np.abs(s.g.xi); hi=xi>0.5*xi.max(); lo=xi<0.25*xi.max()
    msg=f"[{zmin},{zmax}] n={n} dz={s.g.dz:.4f} nsteps={nst} A={A} filt={filt}: "
    for k in range(1,nst+1):
        h=s.step(h)
        if not np.all(np.isfinite(h)): print(msg+f"NaN at s={k*s.ds:.4f}"); return
        if verbose and k%(nst//15)==0:
            sp=np.abs(np.fft.fft(s.wg*h)); G=s.N(h); spG=np.abs(np.fft.fft(s.wg*G))
            # where in z does the high-frequency content of wg*h live?
            hf=np.fft.ifft(np.where(hi,np.fft.fft(s.wg*h),0)); iz=np.argmax(np.abs(hf))
            print(f"   s={k*s.ds:.3f} max|w h|^ low={sp[lo].max():.2e} high={sp[hi].max():.2e}  G^ high={spG[hi].max():.2e}  HF located z={s.g.z[iz]:.1f} amp={np.abs(hf).max():.2e}")
    print(msg+"finite to S; max|G|=%.4e"%np.max(np.abs(s.N(h))))
if __name__=="__main__":
    diag(-64,64,8192,600,0.05,verbose=True)
    for nst in (1200,2400): diag(-64,64,8192,nst,0.05)
    diag(-64,64,8192,600,0.005)
    diag(-64,64,8192,600,0.05,filt=0.5)
    diag(-32,32,4096,600,0.05)
    diag(-32,96,8192,600,0.05)
    diag(-96,32,8192,600,0.05)
