#!/usr/bin/env python3
"""Linchpin: is the tan-Barnes tau ELEMENTARY?  Claim: tau(s)=pi/sin(pi s/2)=Gamma(s/2)Gamma(1-s/2)
solves tau(s+1)/tau(s)=tan(pi s/2), so Ghat_p(s)=Gamma(s+lam)Gamma(s+2)Gamma(s/2)Gamma(1-s/2)
solves the recurrence Ghat(s+1)=R(s)Ghat(s), R(s)=(s+lam)(s+2)tan(pi s/2) -- a product of 4 ordinary
Gammas (Fox-H eigenfunction), with the Wiener-Hopf split Gamma(s/2)[left poles] * Gamma(1-s/2)[right].
If confirmed, the 'unconstructed tan-Barnes' is classical. mpmath, high precision.
"""
import mpmath as mp
mp.mp.dps = 30

def tau(s): return mp.pi / mp.sin(mp.pi * s / 2)
def R(s, lam): return (s + lam) * (s + 2) * mp.tan(mp.pi * s / 2)
def Ghat_p(s, lam): return mp.gamma(s + lam) * mp.gamma(s + 2) * mp.gamma(s / 2) * mp.gamma(1 - s / 2)

print("(1) tau(s+1)/tau(s) =? tan(pi s/2)   [tau=pi/sin(pi s/2)]")
for s in [mp.mpf('0.37'), mp.mpf('1.3')+0.8j, -mp.mpf('0.6')+1.1j, mp.mpf('2.4')]:
    lhs = tau(s + 1) / tau(s); rhs = mp.tan(mp.pi * s / 2)
    print(f"    s={str(s):<16s} ratio={mp.nstr(lhs,8)}  tan={mp.nstr(rhs,8)}  err={mp.nstr(abs(lhs-rhs),3)}")

print("\n(2) Gamma(s/2)Gamma(1-s/2) =? pi/sin(pi s/2)  (reflection => the WH split)")
for s in [mp.mpf('0.37'), mp.mpf('1.3')+0.8j, mp.mpf('2.4')]:
    lhs = mp.gamma(s / 2) * mp.gamma(1 - s / 2); rhs = mp.pi / mp.sin(mp.pi * s / 2)
    print(f"    s={str(s):<16s} GG={mp.nstr(lhs,8)}  pi/sin={mp.nstr(rhs,8)}  err={mp.nstr(abs(lhs-rhs),3)}")

print("\n(3) Ghat_p(s+1)/Ghat_p(s) =? R(s)=(s+lam)(s+2)tan(pi s/2)   [the full recurrence]")
for lam in [mp.mpf('0.7'), mp.mpf('1.5')+1.0j, mp.mpf('2.3')]:
    for s in [mp.mpf('0.41'), mp.mpf('1.2')+0.5j]:
        lhs = Ghat_p(s + 1, lam) / Ghat_p(s, lam); rhs = R(s, lam)
        print(f"    lam={str(lam):<12s} s={str(s):<12s} ratio/R err={mp.nstr(abs(lhs-rhs),3)}")

print("\n(4) decay of Ghat_p on the inversion line Re s=0 (s=i tau): |Ghat_p| ~ exp(-3pi|tau|/2)?")
for lam in [mp.mpf('0.7')]:
    for tau_ in [mp.mpf('2'), mp.mpf('5'), mp.mpf('10')]:
        v = abs(Ghat_p(1j * tau_, lam))
        print(f"    tau={str(tau_):<5s} |Ghat_p|={mp.nstr(v,4)}  exp(-3pi tau/2)={mp.nstr(mp.e**(-3*mp.pi*tau_/2),4)}  ratio={mp.nstr(v/mp.e**(-3*mp.pi*tau_/2),4)}")
print("\n=> if (1)-(3) err~0: the tan-Barnes is ELEMENTARY (4-Gamma product); eigenfunction is Fox-H.")
