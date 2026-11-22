import mpmath as mp, sys
from a5 import Rslit, X
mp.mp.dps=25
print("R_slit applied to X(.;mu) vs closed form 6(c), strip-interior s")
for sig, lam, mu in [(mp.mpf('0.3'), mp.mpc('-0.4','1.3'), mp.mpc('0.5','-0.8')), (mp.mpf('0.4'), mp.mpf('-0.5'), mp.mpf('0.7'))]:
    s = sig+mp.mpc('0.2','0.4')
    r = Rslit(lambda w: X(w, mu), s, lam, sig, T=8)
    c = (X(s, mu)-mp.gamma(1+mu)/mp.gamma(1+lam)*X(s, lam))/(lam-mu)
    print(f"  sigma={sig} lam={mp.nstr(lam,4)} mu={mp.nstr(mu,4)}: rel diff {mp.nstr(abs(r-c)/abs(c),3)}"); sys.stdout.flush()
