import mpmath as mp, sys
mp.mp.dps = 18
N = lambda s: mp.cos(mp.pi*s)*mp.gamma(2*s+2)/mp.pi
def X(s, lam): return mp.gamma(2*s+lam)*mp.gamma(2*s+2)*mp.gamma(s)*mp.gamma(1-s)
def inner(x, s, mu):
    return mp.quad(lambda b: b**mu*(b+x)**(-2*s-mu)/(1+b), [0, x/2, x, 2*x+0.5, 5, 50, mp.inf])
def claimed(s, lam, mu):
    # tau-form with the endpoint singularity removed: Int_0^inf e^{-(lam+2s)tau}[J(1-e^-tau)-J(1)]dtau + J(1)/(lam+2s)
    J1 = inner(mp.mpf(1), s, mu)
    f = lambda tau: mp.e**(-(lam+2*s)*tau)*(inner(1-mp.e**(-tau), s, mu)-J1)
    return N(s)*2*mp.pi*mp.gamma(2*s+mu)*(mp.quad(f, [0, 0.1, 1, 3, 8, 20, 60]) + J1/(lam+2*s))
closed = lambda s, lam, mu: (X(s, mu)-mp.gamma(1+mu)/mp.gamma(1+lam)*X(s, lam))/(lam-mu)
for sig, lam, mu, t in [(0.3, mp.mpc('-0.4','1.3'), mp.mpc('0.5','-0.8'), 0.8), (0.4, mp.mpf('-0.5'), mp.mpf('0.7'), 0.0)]:
    s = mp.mpc(sig, t); a, b = claimed(s, lam, mu), closed(s, lam, mu)
    print(f"6(c) vs claimed (regularised quadrature) sigma={sig} lam={mp.nstr(lam,4)} mu={mp.nstr(mu,4)}: rel diff {mp.nstr(abs(a-b)/abs(b),3)}"); sys.stdout.flush()
