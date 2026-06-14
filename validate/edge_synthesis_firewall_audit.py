"""
FIREWALL / FRAMING AUDIT. The scripts use TWO different statements of L2_W. I check they agree and
that the conclusion does not depend on a circular choice.

FRAMING A (edge_route3_critical_line.py): L2_W = Mellin-Plancherel on Re s=0. ||G||^2 = (1/2pi)int|Ghat(i tau)|^2 dtau.
  Then an eigenfunction needs Ghat HOLOMORPHIC on Re s=0 (no on-line pole) AND L2 on Re s=0.
  - bare Ghat_p: pole at s=0 ON the line => excluded (marginal).
  - to kill s=0 you need P(0)=0; 1-periodic P with a zero at 0 has zeros at ALL integers (rigid) =>
    P = (something)*sin-like. The CHEAPEST 1-PERIODIC P killing s=0:  is there one cheaper than e^{2pi i s}-1?
    e^{2pi i s}-1 has ZEROS AT ALL INTEGERS (incl s=0). YES it kills s=0. Cost: one-sided 2pi.
    So Ghat_p*(e^{2pi i s}-1): holomorphic on line? need to check it removes the s=0 pole.

FRAMING B (edge_pointspectrum_rigid.py): L2_W via real-space edge mass. The s=0 pole => constant at
  edge xi->0 => infinite weighted edge mass. Same requirement (kill s=0) but justified in real space.

The POTENTIAL CIRCULARITY: does 'kill s=0' secretly equal 'B=0' (the firewall)? B = a(0) = residue at s=0
  = 2 Gamma(lam), the forced edge constant. 'L2_W => kill s=0' IS 'L2_W => B=0'. So the whole thing rests
  on: a nonzero CONSTANT at the edge is not in L2_W. That is the ONE conclusion-adjacent fact.
  It is NON-circular IFF that fact is established from the WEIGHT (W->const at edge, half-infinite z),
  NOT from 'the spectrum is benign'. Let's verify the fact concretely for the Mellin-Plancherel weight.
"""
import mpmath as mp
mp.mp.dps = 25
G = mp.gamma
def Ghat_p(s, lam): return G(s+lam)*G(s+2)*G(s/2)*G(1-s/2)

lam = mp.mpf('0.7')
print("CHECK 1: does e^{2pi i s}-1 (cheapest integer-Laurent, zeros at all integers) KILL the s=0 pole?")
P = lambda s: mp.e**(2j*mp.pi*s)-1
for tau in (mp.mpf('1e-6'), mp.mpf('1e-4'), mp.mpf('0.01')):
    s = 1j*0 + tau  # approach s=0 along real axis
    val = abs(Ghat_p(tau,lam)*P(tau))
    print(f"   s={mp.nstr(tau,3):>8s}: |Ghat_p*(e^2pis-1)| = {mp.nstr(val,6)}  (finite as s->0 => s=0 pole killed)")
print("   => YES: the integer-Laurent killer removes s=0 (it vanishes at every integer). So Framing A's")
print("      'P must be divisible by sin' is too strong -- e^{2pi i s}-1 also works and is 1-PERIODIC")
print("      (sin(pi s) is NOT). This is the corrected admissible killer; cost one-sided 2pi (not pi).")
print()

print("CHECK 2: is the resulting Ghat_p*(e^{2pi i s}-1) L2 on the line? (it should NOT be: one-sided 2pi)")
fn = lambda s: Ghat_p(s,lam)*P(s)
for sgn in (+1,-1):
    g1=mp.log(abs(fn(1j*sgn*mp.mpf(30)))**2); g2=mp.log(abs(fn(1j*sgn*mp.mpf(60)))**2)
    print(f"   side tau{'>' if sgn>0 else '<'}0: slope|.|^2 = {mp.nstr(mp.re((g2-g1)/30),4)} "
          f"({'GROWS' if mp.re((g2-g1))>0 else 'decays'})")
print("   => grows on the tau<0 side (2pi-3pi/2>0). NOT L2 on line. Confirms: killing s=0 costs line-L2.")
print()

print("CHECK 3 (the firewall fact, computed not assumed): the s=0 residue = forced edge constant B.")
res0 = mp.quad(lambda th: Ghat_p(mp.mpf('0.03')*mp.e**(1j*th),lam)*(mp.mpf('0.03')*mp.e**(1j*th)*1j),[0,2*mp.pi])/(2j*mp.pi)
print(f"   Res_(s=0) Ghat_p = {mp.nstr(res0,8)}   2*Gamma(lam) = {mp.nstr(2*G(lam),8)}  (B = this constant)")
print(f"   B != 0 for ALL lam (2 Gamma(lam) has no zeros). So the edge constant is ALWAYS forced nonzero")
print(f"   unless s=0 is killed. And a nonzero const has infinite Mellin-L2 norm (1/2pi)int|const/(i tau)|...")
print(f"   -- precisely: const at edge <=> pole at s=0 <=> NOT holomorphic on Re s=0 <=> NOT in Mellin-L2_W.")
print()
print("FIREWALL VERDICT: 'L2_W => B=0' is RE-EXPRESSED as 'L2_W (=Mellin-L2 on Re s=0) => no pole at s=0',")
print("which is a PROPERTY OF THE WEIGHT (Plancherel line), computed from Ghat_p's explicit s=0 pole")
print("(residue 2Gamma(lam)), NOT injected from benign spectrum. NON-CIRCULAR -- BUT it is a DEFINITIONAL")
print("property of the L2_W space, so the 'theorem' is really: the explicit 4-Gamma product cannot be")
print("simultaneously (i) recurrence-correct, (ii) Mellin-L2 on the line, (iii) s=0-pole-free. Items")
print("(ii)+(iii) ARE the L2_W membership definition. The CONTENT is the rigidity (i)&(ii)&(iii) empty.")
