import mpmath as mp, sys
mp.mp.dps = 15
N = lambda s: mp.cos(mp.pi*s)*mp.gamma(2*s+2)/mp.pi
def X(s, lam): return mp.gamma(2*s+lam)*mp.gamma(2*s+2)*mp.gamma(s)*mp.gamma(1-s)
# (ii) 6(c): data X(.;mu) <-> k_mu^ = 2 pi^2 Gamma(2s+mu)/sin(2 pi s) = N^{-1} X ; k_mu = pi * Int_0^inf b^mu u^mu e^{-b u}/(1+b) db
# claimed resolvent: G^(s) = N(s) 2 pi Gamma(2s+mu) Int_0^1 dx (1-x)^{lam+2s-1} Int_0^inf b^mu (b+x)^{-2s-mu}/(1+b) db
print("check N^{-1}X identity:", mp.nstr(abs(X(mp.mpc(0.3,0.7),0.4)/N(mp.mpc(0.3,0.7)) - 2*mp.pi**2*mp.gamma(2*mp.mpc(0.3,0.7)+0.4)/mp.sin(2*mp.pi*mp.mpc(0.3,0.7))),3))
def inner(x, s, mu):
    # Int_0^inf b^mu (b+x)^(-2s-mu)/(1+b) db
    return mp.quad(lambda b: b**mu*(b+x)**(-2*s-mu)/(1+b), [0, x/2, x, 2*x+0.5, 5, 50, mp.inf])
def claimed(s, lam, mu):
    f = lambda x: (1-x)**(lam+2*s-1)*inner(x, s, mu)
    return N(s)*2*mp.pi*mp.gamma(2*s+mu)*mp.quad(f, [0, 0.01, 0.1, 0.5, 0.9, 0.99, 1])
def closed(s, lam, mu):
    return (X(s, mu)-mp.gamma(1+mu)/mp.gamma(1+lam)*X(s, lam))/(lam-mu)
print("(ii) closed form 6(c) vs claimed resolvent on the line Re s = sigma")
for sig, lam, mu in [(0.4, mp.mpf('-0.5'), mp.mpf('0.7')),
                     (0.3, mp.mpc('-0.4','1.3'), mp.mpc('0.5','-0.8')), (0.25, mp.mpc('-0.3','-2.5'), mp.mpf('-0.2')),
                     (0.45, mp.mpc('1.5','6'), mp.mpc('-0.6','0.3'))]:
    worst = 0
    for t in [0.8]:
        s = mp.mpc(sig, t)
        a, b = claimed(s, lam, mu), closed(s, lam, mu)
        worst = max(worst, abs(a-b)/abs(b))
    print(f"  sigma={sig} lam={mp.nstr(lam,4)} mu={mp.nstr(mu,4)}: worst rel diff {mp.nstr(worst,3)}"); sys.stdout.flush()
# (iii) defect: (lam - A) k_lam = pi Gamma(1+lam) * 1 ; check in z-space by quadrature of k_lam and its derivative
def k(u, lam):   # k_lam as function of u=e^{-z/2}
    return mp.pi*mp.quad(lambda b: b**lam*u**lam*mp.e**(-b*u)/(1+b), [0, 1/u, 10/u, mp.inf])
print("(iii) (lam-A)k_lam = lam k + 2 k_z + u k  [k_z = -(u/2) dk/du]  vs  pi*Gamma(1+lam):")
for lam in [mp.mpf('0.3'), mp.mpc('-0.4','1.0')]:
    for u in [mp.mpf('0.2'), mp.mpf(3)]:
        val = lam*k(u, lam) - u*mp.diff(lambda v: k(v, lam), u) + u*k(u, lam)
        print(f"   lam={mp.nstr(lam,3)} u={u}: {mp.nstr(val,12)}   vs {mp.nstr(mp.pi*mp.gamma(1+lam),12)}")
