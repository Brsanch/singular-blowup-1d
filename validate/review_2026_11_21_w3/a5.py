# Attack 5: claimed resolvent  N (lam-A)^{-1} N^{-1}  vs  repo's R_slit (split Cauchy formula, re-implemented from
# validate/edge_slitplane_representation.py) and vs the closed form 6(c).
import mpmath as mp, sys
mp.mp.dps = 30
N = lambda s: mp.cos(mp.pi*s)*mp.gamma(2*s+2)/mp.pi
def X(s, lam): return mp.gamma(2*s+lam)*mp.gamma(2*s+2)*mp.gamma(s)*mp.gamma(1-s)
def Rslit(Fhat, s, lam, sig, T=14):
    Phi = lambda tp: Fhat(sig+1j*tp)/(X(sig+1j*tp, lam)*(lam+2*(sig+1j*tp)))
    wray = lambda tp: mp.e**(4j*mp.pi*sig)*mp.e**(-4*mp.pi*tp)
    def D(w):
        fp = lambda tp: Phi(tp)*(-4*mp.pi*wray(tp))/(wray(tp)-w)
        fm = lambda tp: Phi(tp)*(-4*mp.pi*wray(tp))*(1/(wray(tp)-w)-1/wray(tp))
        return -(mp.quad(fp, [0,0.5,1,2,3,6,T])+mp.quad(fm, [-T,-6,-3,-2,-1,-0.5,0]))/(2j*mp.pi)
    return X(s, lam)*(D(mp.e**(4j*mp.pi*s))-D(mp.mpf(1)))
# claimed: F = N k_b, k_b = u^a e^{-b u}:  Ghat(s) = N(s) 2 Gamma(2s+a) Int_0^1 (1-x)^{lam+2s-1} (b+x)^{-2s-a} dx
def claimed_fam(s, lam, a, b):
    I = mp.quad(lambda x: (1-x)**(lam+2*s-1)*(b+x)**(-2*s-a), mp.linspace(0,1,9))
    return N(s)*2*mp.gamma(2*s+a)*I
if __name__=="__main__": print("(i) family data F^ = N(s) 2 Gamma(2s+a) b^(-2s-a):  claimed resolvent vs R_slit, strip-interior points")
for sig, lam, a, b in [] if __name__!="__main__" else [(mp.mpf('0.4'), mp.mpf('0.3'), 1, mp.mpf(1)),
                       (mp.mpf('0.4'), mp.mpf('-0.5'), 1, mp.mpf(1)),
                       (mp.mpf('0.3'), mp.mpc('-0.35','1.7'), 1, mp.mpc('0.6','0.5')),
                       (mp.mpf('0.2'), mp.mpc('-0.3','-4.0'), mp.mpf('0.5'), mp.mpc('0.2','-0.3')),
                       (mp.mpf('0.45'), mp.mpc('-0.85','0.4'), 2, mp.mpf('3')),
                       (mp.mpf('0.3'), mp.mpc('2.0','10.0'), 1, mp.mpf('0.1'))]:
    Fh = lambda s: N(s)*2*mp.gamma(2*s+a)*b**(-2*s-a)
    worst = 0
    for ds in [mp.mpf('0.07'), mp.mpc('0.2','0.6'), mp.mpc('0.33','-1.1')]:
        s = sig+ds
        r = Rslit(Fh, s, lam, sig); c = claimed_fam(s, lam, a, b)
        worst = max(worst, abs(r-c)/abs(c))
    print(f"  sigma={sig} lam={mp.nstr(lam,4)} a={a} b={mp.nstr(b,3)}: worst rel diff = {mp.nstr(worst,3)}   (|G|~{mp.nstr(abs(c),3)})"); sys.stdout.flush()
