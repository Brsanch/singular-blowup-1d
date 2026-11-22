# strong continuity on a NON-domain datum: g0^ = Gamma(2s+2)/sin(pi s)  (g0^(1/2)=2 != 0), k0 = 1/(1+u)
import mpmath as mp, math
mp.mp.dps=20
N=lambda s: mp.cos(mp.pi*s)*mp.gamma(2*s+2)/mp.pi
def Gh(s,tau):
    if tau==0: return mp.gamma(2*s+2)/mp.sin(mp.pi*s)
    x=1-mp.e**(-tau)
    return N(s)*2*mp.gamma(2*s)*mp.e**(-2*s*tau)*mp.e**x*mp.gammainc(1-2*s,x)
ts=[mp.mpf(j)/8 for j in range(-48,49)]
for sig in (0.3,0.45):
    n0=math.sqrt(sum(abs(Gh(mp.mpc(sig,t),0))**2 for t in ts))
    prev=None
    for tau in (1e-1,1e-2,1e-3,1e-4,1e-5):
        d=math.sqrt(sum(abs(Gh(mp.mpc(sig,t),tau)-Gh(mp.mpc(sig,t),0))**2 for t in ts))/n0
        nt=math.sqrt(sum(abs(Gh(mp.mpc(sig,t),tau))**2 for t in ts))/n0
        print(f"sigma={sig} tau={tau:g}: ||S(tau)g0-g0||/||g0||={d:.5e}  ||S(tau)g0||/||g0||={nt:.6f}"+("" if prev is None else f"  local exponent={math.log10(prev/d):.3f} (1-2sigma={1-2*sig:.2f})"))
        prev=d
