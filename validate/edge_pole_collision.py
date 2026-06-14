#!/usr/bin/env python3
"""
ROUTE 1 — residue / pole-collision closure of the BALANCE step (§5.4 of edge_mellin_symbol.md).
Sharpens the heuristic "killing both bad classes costs e^{2pi|tau|}" into a structural lower bound.

Eigenfunction Mellin transform: Ghat_p(s) = G(s+lam)G(s+2)G(s/2)G(1-s/2); general solution
Ghat = Ghat_p * P with P 1-periodic. L2_W eigenfunction <=> G(xi) decays at BOTH ends <=>
Ghat_p*P holomorphic on the strip between the edge & far-field indicial lines AND L2 on Re s=0.

The argument (all spectrum-free, firewall-clean):
  (1) EDGE residue: Res_{s=2} Ghat_p = -2 G(2+lam)G(4)G(1) = -12 Gamma(2+lam), NONZERO for all
      generic lam => with P=1 the edge term xi^{-2} diverges (xi->0). (analytic + numeric).
  (2) POLE ORDERS: right poles s=2,4,.. SIMPLE (class {0} mod 1); left poles s=0 simple (marginal,
      in N), s=-2,-4 DOUBLE (G(s/2)&G(s+2)) (class {0}), s=-lam-k SIMPLE (class {-lam mod 1}).
  (3) STRUCTURAL THEOREM (closes the meromorphic-P loophole): a 1-periodic meromorphic P has ALL
      its poles/zeros with Re s in (-1/2,1/2] (one period). Ghat_p is ZERO-FREE (product of Gammas),
      so any pole of P lands on/beside the inversion line and SURVIVES in Ghat_p*P => not L2.
      Hence P must be ENTIRE: P(s)=sum_n c_n e^{2pi i n s} (Laurent in w=e^{2pi i s}).
  (4) GROWTH PER ZERO: an entire 1-periodic P with D zeros per period (mult.) has minimal symmetric
      two-sided line-growth e^{pi D |tau|} (centered Fourier support n in [-D/2,D/2]). [sin: D=1 -> e^{pi|tau|}.]
  (5) REQUIRED zeros to make Ghat_p*P pole-free on the strip (decay both ends): integer class {0}
      needs order max(edge 1, far-double 2)=2; class {-lam} needs order 1. D_req >= 3 (generic lam).
      [Even the over-optimistic D=2 ignoring the far double pole already loses.]
  (6) INEQUALITY: line-slope(Ghat_p*P) = -3pi/2 + pi*D_req >= -3pi/2 + 2pi = +pi/2 > 0  (D>=2),
      or +3pi/2 (D=3). POSITIVE => NOT L2 => NO eigenfunction, for every non-integer Re lam>=0.
  (7) COLLISION: lam in Z merges {-lam} into {0}; those modes are lam=1 (symmetry, in N) and lam
      even (tan=0 => REGULAR, lam<=0). No Re lam>0 point eigenvalue mod N.

REMAINING GAP (honestly stated): step (4) gives the minimal growth for ENTIRE 1-periodic P whose
zeros are exactly the required classes; a full Pincherle proof must also exclude entire-P with
EXTRA zeros arranged to lower one-sided growth while still controlling the other end (the centering
argument bounds this but a clean theorem that no entire 1-periodic factor beats e^{pi D|tau|} for
TWO-SIDED L2 is the last step). The structural theorem (3) + the D>=2 margin make the conclusion
robust, but "rigorous" is reserved until (4) is a stated lemma with proof.
"""
import mpmath as mp
mp.mp.dps = 30
G = mp.gamma

def Ghat_p(s, lam): return G(s+lam)*G(s+2)*G(s/2)*G(1-s/2)
def R(s, lam):      return (s+lam)*(s+2)*mp.tan(mp.pi*s/2)

def res_at(s0, lam, eps=mp.mpf('1e-5')):
    f = lambda th: Ghat_p(s0+eps*mp.e**(1j*th), lam)*(eps*mp.e**(1j*th)*1j)
    return mp.quad(f, [0, 2*mp.pi])/(2j*mp.pi)

def pole_order(s0, lam):
    for k in range(0, 5):
        v0 = (mp.mpf('1e-4')**k)*abs(Ghat_p(s0+mp.mpf('1e-4'), lam))
        v1 = (mp.mpf('1e-6')**k)*abs(Ghat_p(s0+mp.mpf('1e-6'), lam))
        if v0 > 1e-25 and abs(mp.log(v1/v0)) < 1:
            return k
    return -1

def line_slope(f, t1=mp.mpf(30), t2=mp.mpf(60)):
    return (mp.log(abs(f(1j*t2)))-mp.log(abs(f(1j*t1))))/(t2-t1)

print("(1) EDGE residue Res_{s=2} Ghat_p (xi^-2 term):  -12*Gamma(2+lam), nonzero for all lam")
for lam in [mp.mpf('0.7'), mp.mpf('2.5'), mp.mpf('0.5')]:
    num = res_at(2, lam); ana = -12*G(2+lam)
    print(f"    lam={str(lam):<5s} numeric={mp.nstr(num,8)}  analytic -12*Gamma(2+lam)={mp.nstr(ana,8)}"
          f"  err={mp.nstr(abs(num-ana),3)}")

print("\n(2) pole orders of Ghat_p  (right simple class{0}; left s=-2,-4 DOUBLE class{0}; s=-lam SIMPLE):")
for lam in [mp.mpf('0.7'), mp.mpf('1')]:
    print(f"    lam={lam}: s=2:{pole_order(mp.mpf(2),lam)}  s=0:{pole_order(mp.mpf(0),lam)} "
          f" s=-2:{pole_order(mp.mpf(-2),lam)}  s=-lam:{pole_order(-lam,lam)}")

print("\n(3) STRUCTURAL: a real-positive pole of P (w=wp>0) sits at Re s=0 (ON contour) -> P must be ENTIRE.")
wp = mp.mpf('1e6'); s_p_re = mp.arg(wp)/(2*mp.pi); s_p_im = -mp.log(abs(wp))/(2*mp.pi)
print(f"    pole at w={mp.nstr(wp,2)} -> s={mp.nstr(s_p_re,3)}+{mp.nstr(s_p_im,3)}i  (Re s=0, on the line)")

print("\n(4) per-zero growth: sin(pi s) has D=1 zero/period, line-slope:")
print(f"    measured {mp.nstr(line_slope(lambda s: mp.sin(mp.pi*s)),6)}  target +pi={mp.nstr(mp.pi,6)}")

print("\n(5)+(6) minimal-killer product slope (D=2 over-optimistic, D=3 correct); want <0 for L2:")
print(f"    Ghat_p decay slope = {mp.nstr(line_slope(lambda s: Ghat_p(s, mp.mpf('0.7'))),6)} (-> -3pi/2)")
print(f"    {'lam':>10s}  {'D=2 slope':>10s}  {'D=3 slope':>10s}  verdict")
for lam in [mp.mpf('0.05'), mp.mpf('0.3'), mp.mpf('0.7'), mp.mpf('1.4'), mp.mpf('2.5'),
            mp.mpf('3.6'), mp.mpf('0.5')+1.5j]:
    fD2 = lambda s: Ghat_p(s, lam)*mp.sin(mp.pi*s)*mp.sin(mp.pi*(s+lam))
    fD3 = lambda s: Ghat_p(s, lam)*mp.sin(mp.pi*s)**2*mp.sin(mp.pi*(s+lam))
    s2 = line_slope(fD2); s3 = line_slope(fD3)
    pos = (s2.real if hasattr(s2, 'real') else s2) > 0
    print(f"    {str(lam):>10s}  {mp.nstr(s2,5):>10s}  {mp.nstr(s3,5):>10s}  {'NOT L2' if pos else '??'}")

print("\n(A) one-sided (D=1) killer is INSUFFICIENT: Ghat_p*sin(pi s) is L2-on-line but the {-lam}")
print("    far-field pole survives (G(xi)~xi^lam at xi->inf):")
lam = mp.mpf('0.7')
for e in [mp.mpf('1e-3'), mp.mpf('1e-5'), mp.mpf('1e-7')]:
    v = abs(Ghat_p(-lam+e, lam)*mp.sin(mp.pi*(-lam+e)))
    print(f"    eps={mp.nstr(e,2)}: |Ghat_p*sin| near s=-lam = {mp.nstr(v,6)}  (~1/eps => pole survives)")

print("\n=> ROUTE 1: no Re lam>=0 point eigenvalue mod N (D_req*pi >= 2pi > 3pi/2). Firewall: the two")
print("   exponents (-3pi/2 Stirling, +pi/zero periodicity) are spectrum-free; the inequality is arithmetic.")
print("   Remaining gap: a clean Pincherle lemma that no entire 1-periodic P beats e^{pi D|tau|} two-sided.")
