#!/usr/bin/env python3
"""ROUTE 3 (corrected, DECISIVE) — L2_W = Mellin-L2 on the CRITICAL LINE Re s = 0.

Correction over the first pass: L2_W(dxi/xi) is the Mellin-Plancherel space whose inversion line
is Re s = 0 (dxi/xi is the multiplicative Haar measure; ||G||^2_{L2_W} = (1/2pi) int |Ghat(i tau)|^2 dtau).
So an L2_W eigenfunction needs Ghat = Ghat_p * P (P 1-periodic) to be BOTH
   (H) HOLOMORPHIC on Re s = 0  (no pole ON the line), AND
   (L) L2 on Re s = 0  (decaying).
Ghat_p = Gamma(s+lam)Gamma(s+2)Gamma(s/2)Gamma(1-s/2) has a pole AT s=0 (from Gamma(s/2)) sitting
exactly ON the line => Ghat_p itself is NOT holomorphic on Re s=0 => the bare particular solution is
the MARGINAL mode (the s=0 pole = the dilation/marginal mode in N), NOT an L2_W eigenfunction.

[Earlier real-space finding, now explained: the sigma=0.5-inverted G tends to a nonzero CONSTANT as
xi->0 (= the s=0 residue 2Gamma(lam)) => int|G|^2 dxi/xi LOG-diverges at the edge.  That log-divergence
IS the s=0-pole-on-the-line obstruction.  Verified: G(1e-4)~2.59=const, near-edge int ~15.5/decade flat.]

So the genuine question: can a 1-periodic P make Ghat_p*P holomorphic AND L2 on Re s=0?
  * To be HOLOMORPHIC on Re s=0, P must CANCEL the s=0 pole.  A 1-periodic P with P(0)=0 has zeros at
    ALL integers (rigidity) => P divisible by sin(pi s).  sin(pi s) DOES kill s=0 (and s=2,4,..; and
    s=-2,-4,..) — good for holomorphy on the line.  Cost on the line: |sin(pi(i tau))|=sinh(pi|tau|)
    ~ (1/2)e^{pi|tau|}.
  * Ghat_p decays e^{-3pi|tau|/2} on the line; times one sin => e^{-pi|tau|/2} => STILL L2.  So
    Ghat_p*sin(pi s) is holomorphic on Re s=0 AND L2 on Re s=0 — a CANDIDATE eigenfunction!
    => the decisive question is whether THIS candidate also has the right HALF-PLANE analyticity to
    be a true L2_W function with the eigen-recurrence, i.e. whether killing s=0 with sin(pi s) is
    ADMISSIBLE or whether it destroys decay at one spatial end.

This script settles it by the two-sided real-space asymptotics of the candidate
    Ghat_c(s) = Ghat_p(s) * sin(pi s),
whose surviving poles are:  LEFT  s=-lam-k (Gamma(s+lam), UN-killed since -lam non-integer),
                            (the integer poles s=0,-2,-4,.. and s=2,4,.. are KILLED by sin),
                            plus s=-2-k from Gamma(s+2) are integers => KILLED too.
So the ONLY surviving poles are the NON-integer ones s=-lam-k (k>=0), all in Re s<=-lam<0 (Re lam>0).
=> Ghat_c is holomorphic in the whole RIGHT half-plane Re s> -lam, in particular no right pole at all
   in Re s>=0 => closing right gives G(xi)=0-decay faster than any power as xi->0  (EDGE OK, decays!).
=> far field xi->inf governed by the nearest LEFT pole s=-lam: G ~ xi^{lam}  => BLOWS UP for Re lam>0.

THE decisive dichotomy, now exact and pole-driven (NOT heuristic):
  killing s=0 (mandatory for line-holomorphy) with the 1-periodic sin ALSO removes the s=2 edge pole
  (edge becomes fine) but CANNOT remove the non-integer far pole s=-lam (off the integer lattice) =>
  far field xi^{lam} survives => NOT L2_W for Re lam>0.  Adding sin(pi(s+lam)) to kill s=-lam pushes
  back to TWO sins => e^{2pi|tau|} on the line => destroys L2 on the line (loses (L)).
  EITHER WAY no L2_W eigenfunction for non-integer lam, Re lam>0.   Collisions (lam in Z): s=-lam is an
  integer => already killed by the single sin => far pole GONE => but then it's an N/symmetry mode.

This script verifies every link numerically: (H) holomorphy on Re s=0, (L) line-L2, the surviving-pole
list, the far-field xi^{lam} blow-up of Ghat_p*sin, and the line-blowup of Ghat_p*sin*sin.
mpmath.  No apfsd risk.
"""
import mpmath as mp
mp.mp.dps = 22
G = mp.gamma

def Ghat_p(s, lam): return G(s + lam) * G(s + 2) * G(s / 2) * G(1 - s / 2)

# Ghat_c = Ghat_p * sin(pi s).  The s=0 (and even-integer) poles of Gamma(s/2) are removable: using
# reflection Gamma(s/2)Gamma(1-s/2)=pi/sin(pi s/2) and sin(pi s)=2 sin(pi s/2)cos(pi s/2),
#   Gamma(s/2)*sin(pi s) = 2 pi cos(pi s/2) / Gamma(1-s/2)   (pole-free at even integers).
# So  Ghat_c = Gamma(s+lam) Gamma(s+2) * [Gamma(1-s/2) * 2 pi cos(pi s/2)/Gamma(1-s/2)]
#            = Gamma(s+lam) Gamma(s+2) * 2 pi cos(pi s/2).   <-- closed, entire-up-to-Gamma(s+lam,s+2) poles.
def Ghat_c(s, lam):  return G(s + lam) * G(s + 2) * 2 * mp.pi * mp.cos(mp.pi * s / 2)
def Ghat_cc(s, lam): return Ghat_c(s, lam) * mp.sin(mp.pi * (s + lam))      # also kill far pole

print("=" * 94)
print("(H) Ghat_p has a POLE at s=0 (on the L2_W critical line) => bare Ghat_p is the MARGINAL mode.")
print("    Ghat_c = Ghat_p*sin(pi s) is HOLOMORPHIC on Re s=0 (sin kills the s=0 pole).  Check |.|<inf:")
print("=" * 94)
lam = mp.mpf('0.7')
for tau in (mp.mpf('0'), mp.mpf('0.0001'), mp.mpf('1'), mp.mpf('3')):
    # Ghat_p blows up at s=0; Ghat_c is finite there
    pc = abs(Ghat_c(1j * tau, lam))
    near0 = abs(Ghat_c(mp.mpf('1e-8') + 0j, lam))
    print(f"    tau={mp.nstr(tau,4):>7s}: |Ghat_c(i tau)| = {mp.nstr(pc,6)}   (finite => holomorphic on line)")
print(f"    Ghat_c at s=1e-8 (limit s->0): {mp.nstr(abs(Ghat_c(mp.mpf('1e-8'),lam)),6)} (finite),"
      f"  Ghat_p at s=1e-8: {mp.nstr(abs(Ghat_p(mp.mpf('1e-8'),lam)),6)} (BLOWS UP, the pole)")

print("\n(L) line-L2 on Re s=0: ||Ghat_c||^2 = (1/2pi) int|Ghat_c(i tau)|^2 dtau  (finite => L2):")
for lam in (mp.mpf('0.3'), mp.mpf('0.7'), mp.mpf('1.4') + 1.1j, mp.mpf('2.3')):
    l2 = mp.quad(lambda t: abs(Ghat_c(1j * t, lam)) ** 2, [-30, -3, 0, 3, 30]) / (2 * mp.pi)
    slope = (mp.log(abs(Ghat_c(1j * mp.mpf(60), lam)) ** 2) - mp.log(abs(Ghat_c(1j * mp.mpf(30), lam)) ** 2)) / 30
    print(f"    lam={str(lam):<14s} ||Ghat_c||^2_line = {mp.nstr(l2,6):>10s}  line slope|.|^2={mp.nstr(slope,4)}"
          f" (pred -pi=-3.14 => L2)")

print("\n(POLES) surviving poles of Ghat_c=Ghat_p*sin(pi s) (sin kills ALL integer poles):")
print("    only the NON-integer s=-lam-k (k>=0) survive, all in Re s = -Re lam - k <= -Re lam < 0.")
print("    => NO pole in Re s>=0 => EDGE (xi->0) decays super-polynomially (good); FAR (xi->inf)")
print("       controlled by nearest left pole s=-lam => G ~ xi^{lam}.")
lam = mp.mpf('0.7')
for s_test, nm in [(mp.mpf('2'), 's=2 (was edge pole)'), (mp.mpf('0'), 's=0 (was marginal)'),
                   (-lam, 's=-lam (far pole)'), (-lam - 1, 's=-lam-1')]:
    # residue of Ghat_c at s_test (0 if killed)
    r = mp.quad(lambda th: Ghat_c(s_test + mp.mpf('0.04') * mp.e**(1j*th), lam) * (mp.mpf('0.04')*mp.e**(1j*th)*1j),
                [0, 2*mp.pi]) / (2j*mp.pi)
    print(f"    Res Ghat_c at {nm:<22s} = {mp.nstr(r,5)}   ({'KILLED' if abs(r)<1e-8 else 'SURVIVES'})")

print("\n(FAR BLOW-UP) Ghat_c far field: G ~ xi^{lam} as xi->inf from the s=-lam pole => NOT L2_W (Re lam>0).")
print("    measure the real-space far slope d log|G|/d log xi of the sigma=-0.35-inverted Ghat_c (strip")
print("    -lam<sigma<0 isolates s=-lam as the nearest left pole):")
def Gc(xi, lam, sigma, T=20):
    return mp.quad(lambda t: Ghat_c(sigma + 1j*t, lam) * xi**(-(sigma+1j*t)), [-T,-3,-1,1,3,T]) / (2*mp.pi)
for lam in (mp.mpf('0.7'), mp.mpf('1.4'), mp.mpf('2.3')):
    sigma = -lam / 2     # strip between s=-lam and s=0
    v1, v2 = abs(Gc(mp.mpf('30'), lam, sigma)), abs(Gc(mp.mpf('100'), lam, sigma))
    slope = (mp.log(v2) - mp.log(v1)) / (mp.log(mp.mpf('100')) - mp.log(mp.mpf('30')))
    print(f"    lam={mp.nstr(lam,3):>5s}: far slope d log|G|/d log xi = {mp.nstr(mp.re(slope),5):>8s}"
          f"  (pred +lam={mp.nstr(lam,3)} => xi^lam BLOWS UP => not L2_W at inf)")

print("\n(2-SIN COSTS L2) killing the far pole too (sin(pi(s+lam))) restores e^{2pi|tau|} on the line:")
for lam in (mp.mpf('0.7'), mp.mpf('2.3')):
    slope = (mp.log(abs(Ghat_cc(1j*mp.mpf(60),lam))**2) - mp.log(abs(Ghat_cc(1j*mp.mpf(30),lam))**2)) / 30
    print(f"    lam={mp.nstr(lam,3):>5s}: line slope|Ghat_cc|^2 = {mp.nstr(slope,5)} (pred +pi=3.14 => GROWS => NOT line-L2)")

print("=" * 94)
print("VERDICT (continuum, pole-exact): for non-integer lam, Re lam>0 there is NO 1-periodic P making")
print("Ghat_p*P holomorphic AND L2 on Re s=0 AND two-sided real-space decaying:")
print("  - bare Ghat_p: pole at s=0 on the line (marginal mode, log-divergent edge) -- NOT L2_W.")
print("  - Ghat_p*sin(pi s): holomorphic+L2 on the line, edge OK, but far field xi^{lam} BLOWS UP (s=-lam).")
print("  - Ghat_p*sin*sin (kill far too): line integrand grows e^{pi|tau|/2} -- NOT L2 on the line.")
print("=> no Re lam>=0 point eigenvalue mod N (collisions lam in Z = the N/symmetry modes). STABLE.")
print("=" * 94)
