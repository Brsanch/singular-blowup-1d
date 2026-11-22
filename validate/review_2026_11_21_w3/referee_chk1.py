import numpy as np, mpmath as mp
from scipy.special import loggamma
mp.mp.dps=30
N=lambda s: mp.cos(mp.pi*s)*mp.gamma(2*s+2)/mp.pi
c=lambda w: -mp.tan(mp.pi*w)/(2*w+1)
# 1 symbol identity
err=0
for s in [mp.mpc(0.3,0.7),mp.mpc(-0.2,3),mp.mpc(0.1,-11),mp.mpc(-0.45,0.01)]:
    err=max(err,abs(N(s+0.5)*c(s+0.5)+N(s))/abs(N(s)))
print("symbol identity rel err",err)
# two-sided bound: ratio |N(x+it)|/(1+|t|)^(2x+3/2), log-space
def logabsN(x,t):
    s=x+1j*t
    # log|cos(pi s)| = 0.5*log(cos^2(pi x)+sinh^2(pi t))
    lc=np.where(np.abs(t)<20,0.5*np.log(np.cos(np.pi*x)**2+np.sinh(np.pi*np.clip(t,-20,20))**2),np.pi*np.abs(t)-np.log(2))
    return lc+np.real(loggamma(2*s+2))-np.log(np.pi)
t=np.concatenate([np.linspace(-50,50,20001),np.logspace(2,6,400),-np.logspace(2,6,400)])
for x in [-0.49,-0.4,-0.2,0,0.2,0.4,0.49]:
    r=np.exp(logabsN(x,t)-(2*x+1.5)*np.log1p(np.abs(t)))
    lim=np.sqrt(2*np.pi)*2**(2*x+1.5)/(2*np.pi)
    print(f"x={x:+.2f} min ratio {r.min():.4f} max {r.max():.4f}  Stirling limit {lim:.4f} ratio@1e6 {r[np.argmax(t)]:.4f}")
# 3 sign of Re c
for x in [0.05,0.2,0.3,0.45,0.499,0.5]:
    tt=np.linspace(-30,30,60001); tt=tt[np.abs(tt)>1e-9]
    w=x+1j*tt; v=np.real(-np.tan(np.pi*w)/(2*w+1))
    print(f"Re c on Re w={x}: max {v.max():.3e}")
for x in [0.51,0.6,0.75,0.9]:
    tt=np.linspace(-3,3,60001); w=x+1j*tt; v=np.real(-np.tan(np.pi*w)/(2*w+1))
    pos=tt[v>0]
    print(f"Re c on Re w={x}: max {v.max():.3e}, positive for |t|<{np.abs(pos).max():.4f}")
# 2 defect constant: k_lam(z)= (1/2pi) int e^{(sig+it)z} khat dt ; check lam k +2k' + u k = pi Gamma(1+lam)
def chk(sig,lam):
    kh=lambda s: 2*mp.pi**2*mp.gamma(2*s+lam)/mp.sin(2*mp.pi*s)
    def k(z,d=0):
        f=lambda t: mp.exp((sig+1j*t)*z)*kh(sig+1j*t)*(sig+1j*t)**d
        return mp.quad(f,mp.linspace(-12,12,25))/(2*mp.pi)
    for z in [-3,-0.5,1.0,4.0]:
        val=lam*k(z)+2*k(z,1)+mp.exp(-mp.mpf(z)/2)*k(z)
        print(f" sig={sig} lam={lam} z={z}: (lam-A)k={mp.nstr(val,12)}  piGamma(1+lam)={mp.nstr(mp.pi*mp.gamma(1+lam),12)}")
chk(mp.mpf('0.3'),mp.mpc(0.4,0.0)); chk(mp.mpf('0.2'),mp.mpc(-0.1,1.3))
