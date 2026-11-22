# Attack 1: identity, zero-freeness, two-sided symbol bounds
import numpy as np, mpmath as mp
from scipy.special import loggamma
def logN(s):
    t=s.imag
    lc=np.where(t>=0, -1j*np.pi*s-np.log(2)+np.log1p(np.exp(2j*np.pi*s*(t>=0))), 1j*np.pi*s-np.log(2)+np.log1p(np.exp(-2j*np.pi*s*(t<0))))
    return lc+loggamma(2*s+2)-np.log(np.pi)
# identity with mpmath
mp.mp.dps=30
N=lambda s: mp.cos(mp.pi*s)*mp.gamma(2*s+2)/mp.pi
c=lambda s: -mp.tan(mp.pi*s)/(2*s+1)
worst=0
for s in [mp.mpc(0.1,0.3),mp.mpc(-0.4,7),mp.mpc(0.45,-20),mp.mpc(-0.2,0.01),mp.mpc(0.3,50)]:
    worst=max(worst,abs(N(s+0.5)*c(s+0.5)+N(s))/abs(N(s)))
print("identity worst rel",mp.nstr(worst,3))
t=np.concatenate([np.linspace(-1e4,1e4,400001),np.linspace(-5,5,100001)])
print("sigma'   min|N|/(1+|t|)^e   max|N|/(1+|t|)^e   cond   |N(sig')| at t=0   fitted exponent(t=1e3..1e4)")
for sp in [-0.49,-0.45,-0.3,-0.2,-0.1,-0.01,0.01,0.1,0.2,0.3,0.4,0.45,0.49]:
    s=sp+1j*t; e=2*sp+1.5
    r=np.exp(logN(s).real-e*np.log1p(np.abs(t)))
    l1=logN(np.array([sp+1e3j])).real[0]; l2=logN(np.array([sp+1e4j])).real[0]
    print(f"{sp:+.2f}  {r.min():.4e}  {r.max():.4e}  {r.max()/r.min():.3e}  {abs(complex(N(sp))):.4e}  {(l2-l1)/np.log(10):.5f} (claim {e:.2f})")
# zero-free check in the open strip: argument principle on rectangle [-0.499,0.499]x[-30,30]
f=lambda s: mp.diff(N,s)/N(s)
pts=[mp.mpc(-0.499,-30),mp.mpc(0.499,-30),mp.mpc(0.499,30),mp.mpc(-0.499,30),mp.mpc(-0.499,-30)]
mp.mp.dps=20
tot=0
for a,b in zip(pts[:-1],pts[1:]):
    tot+=mp.quad(lambda x: f(a+(b-a)*x)*(b-a),mp.linspace(0,1,61))
print("zeros-poles in rectangle:",mp.nstr(tot/(2j*mp.pi),6))
