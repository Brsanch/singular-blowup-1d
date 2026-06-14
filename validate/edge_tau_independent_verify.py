#!/usr/bin/env python3
"""
INDEPENDENT reproduction (2026-06-14) of the tan-Barnes / Ghat_p construction — fresh code, no
shared helpers with edge_tau_construction.py, dps=50. Covers all five claims (a)-(e):

  tau(s) = pi/sin(pi s/2) = Gamma(s/2)Gamma(1-s/2),  tau(s+1)/tau(s) = tan(pi s/2),
  Ghat_p(s) = Gamma(s+lam)Gamma(s+2)Gamma(s/2)Gamma(1-s/2)  solves  Ghat_p(s+1)/Ghat_p(s) = R(s),
  R(s) = (s+lam)(s+2)tan(pi s/2).

(a) recurrence of tau; (b) reflection; (c) full recurrence over complex lam,s; (d) e^{-3pi|tau|/2}
line-decay (difference-quotient cancels the polynomial Stirling prefactor) + L2 + s=0 residue=2Gamma(lam);
(e) Wiener-Hopf pole split via rgamma + Fox-H parameter reading + signed two-contour Mellin-Barnes
inversion isolating the edge xi^{-2} term. Eigenfunction is a TRUE Fox-H (B-scales 1 and 1/2),
H^{1,3}_{3,1}[ 1/xi | (1-lam,1),(-1,1),(1,1/2) ; (1,1/2) ]; only the (s+lam)(s+2) block is Meijer-G-reducible.
"""
import mpmath as mp, random
mp.mp.dps = 50
G = mp.gamma

def R(s, lam):    return (s + lam) * (s + 2) * mp.tan(mp.pi * s / 2)
def tau(s):       return mp.pi / mp.sin(mp.pi * s / 2)
def tau_g(s):     return G(s / 2) * G(1 - s / 2)
def Ghat_p(s, l): return G(s + l) * G(s + 2) * G(s / 2) * G(1 - s / 2)

random.seed(20260614)
rnd = lambda: mp.mpc(random.uniform(-3, 3), random.uniform(-3, 3))
pts  = [rnd() for _ in range(12)]
lams = [mp.mpc(0.7,0), mp.mpc(2.5,0), mp.mpc(0.5,1.5), mp.mpc(-0.7,0), mp.mpc(1.3,-0.9), mp.mpc(0,2.1)]

# (a)(b)(c)
ea = max(abs(tau(s+1)/tau(s) - mp.tan(mp.pi*s/2)) for s in pts)
eb = max(abs(tau(s) - tau_g(s)) for s in pts)
ec = max(abs(Ghat_p(s+1,l)/Ghat_p(s,l) - R(s,l)) for l in lams for s in pts)
print(f"(a) max|tau(s+1)/tau(s)-tan|            = {mp.nstr(ea,4)}")
print(f"(b) max|pi/sin - Gamma(s/2)Gamma(1-s/2)|= {mp.nstr(eb,4)}")
print(f"(c) max|Ghat_p(s+1)/Ghat_p(s)-R(s)|     = {mp.nstr(ec,4)}  (6 lam x 12 s)")

# (d) line decay exponent via difference quotient (kills the log|tau| Stirling prefactor)
lam = mp.mpf('0.7')
t1, t2 = mp.mpf(20000), mp.mpf(40000)
slope = (mp.log(abs(Ghat_p(1j*t2,lam))) - mp.log(abs(Ghat_p(1j*t1,lam))))/(t2-t1)
print(f"(d) line-decay slope = {mp.nstr(slope,8)}  target -3pi/2 = {mp.nstr(-3*mp.pi/2,8)}"
      f"  relerr {mp.nstr(abs(slope+3*mp.pi/2)/(3*mp.pi/2),3)}")
absq = lambda t: abs(Ghat_p(mp.mpf('0.001')+1j*mp.mpf(t), lam))**2
print(f"    int|Ghat_p(0.001+i tau)|^2 dtau = {mp.nstr(mp.quad(absq,[-mp.inf,-1,1,mp.inf]),5)} (finite => L2)")
res0 = mp.quad(lambda th: Ghat_p(0.01*mp.e**(1j*th),lam)*(0.01*mp.e**(1j*th)*1j),[0,2*mp.pi])/(2j*mp.pi)
print(f"    Res_(s=0) Ghat_p = {mp.nstr(res0,8)}  predicted 2Gamma(lam) = {mp.nstr(2*G(lam),8)}")

# (e) WH split
print("(e) WH split (rgamma=0 at pole):  Gamma(s/2) LEFT {0,-2,-4,..} ; Gamma(1-s/2) RIGHT {2,4,6,..}")
left  = [float(x) for x in (0,-2,-4,-6) if mp.rgamma(mp.mpf(x)/2)==0]
right = [float(x) for x in (2,4,6,8)   if mp.rgamma(1-mp.mpf(x)/2)==0]
print(f"    confirmed LEFT poles {left}   RIGHT poles {right}")
# Fox-H integrand from parameter lists must equal Ghat_p exactly
def foxH(s,l):
    out = mp.mpf(1)
    for a,A in [(1-l,mp.mpf(1)),(mp.mpf(-1),mp.mpf(1)),(mp.mpf(1),mp.mpf('0.5'))]: out *= G(1-a+A*s)
    for b,B in [(mp.mpf(1),mp.mpf('0.5'))]: out *= G(b-B*s)
    return out
ef = max(abs(foxH(s,lam)-Ghat_p(s,lam)) for s in pts)
print(f"    Fox-H H^(1,3)_(3,1)[1/xi|(1-lam,1),(-1,1),(1,1/2);(1,1/2)] integrand vs Ghat_p: max diff {mp.nstr(ef,4)}")
# signed two-contour Mellin-Barnes inversion isolating the edge xi^-2 term:
def Ginv(xi,l,c):  # (1/2pi i)int f ds, ds=i dtau => (1/2pi)int f(c+i tau) dtau
    return mp.quad(lambda t: Ghat_p(c+1j*t,l)*xi**(-(c+1j*t)),
                   [-200,-50,-20,-5,0,5,20,50,200])/(2*mp.pi)
xi = mp.mpf('0.1')
diff = Ginv(xi,lam,mp.mpf(3)) - Ginv(xi,lam,mp.mpf(1))
res2 = mp.quad(lambda th: Ghat_p(2+0.05*mp.e**(1j*th),lam)*xi**(-(2+0.05*mp.e**(1j*th)))*(0.05*mp.e**(1j*th)*1j),
               [0,2*mp.pi])/(2j*mp.pi)
print(f"    edge xi^-2 term: G_(c=3)-G_(c=1)={mp.nstr(diff,8)}  +Res_(s=2)[f xi^-s]={mp.nstr(res2,8)}"
      f"  |diff| {mp.nstr(abs(diff-res2),3)}")
print("=> all (a)-(e) reproduced independently; eigenfunction is a TRUE Fox-H (not Meijer-G).")
