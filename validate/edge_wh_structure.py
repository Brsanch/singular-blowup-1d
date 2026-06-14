#!/usr/bin/env python3
"""Mellin-line structure of the Conjecture-2.4 eigenvalue problem (continuum; immune to the
discrete-Hilbert bias). Eigenproblem: Ghat(s+1)=R(s)Ghat(s), R(s)=(s+lambda)(s+2)tan(pi s/2).

Two probes (verified workflow wf_2d8a4427):
  (D) LATTICE-DIVERGENCE: the terminating-up power series G(xi)=sum a_m xi^m (top power M=-lambda,
      infinite down-tail) blows up at the edge xi->0 like exp(c/sqrt xi) with c ~lambda-INDEPENDENT,
      for every lambda (incl. complex, negative). => no lattice/power-series solution is L^2_W off
      Re=0; the genuine eigenfunctions are MELLIN-LINE integrals (termination is the wrong object).
  (W) WINDING along Re s=0: W(lambda)=(1/2pi)[arg R(i tau)]_{-inf}^{+inf}. Decomposes EXACTLY per
      factor: W = sign(Re lambda)/2 + 1/2 (s+2) + 1/2 (tan). So W=+3/2 (Re lambda>0, constant) and
      +1/2 (Re lambda<0): jump = EXACTLY +1 at Re lambda=0 (one s=-lambda zero crossing the line).

CRITICAL HONESTY (this is NOT a no-eigenvalue certificate):
  * R(0)=0 on Re s=0 for EVERY lambda (tan(0)=0) => the operator is NON-FREDHOLM on the Plancherel
    line (Re s=0 is in the essential spectrum); a Fredholm index is not even well-defined there, and
    the winding is contour-side-ambiguous by 1 (right of s=0: 3/2; left: 1/2).
  * W=3/2 is the RAW symbol winding, NOT the Fredholm index (the (s+lambda)(s+2) piece is the
    entire-function Gamma-growth factor, not a bounded Toeplitz symbol).
  * Even a well-defined index 0 does NOT give dim ker=0: Coburn's lemma is the upgrade and it is
    BLOCKED by the on-line zero. W(lambda) is essentially the indicial-root tautology
    Re(-lambda)<0 <=> Re lambda>0 -- genuinely spectrum-independent (firewall clean) but for the
    SAME reason carries no no-eigenvalue content. The genuine spectral content is the tan-Barnes
    tau(s) (tau(s+1)/tau(s)=tan(pi s/2)), holomorphic & L^2-decaying on the contour -- UNCONSTRUCTED.
Status: reduced-to-sharp-WH-factorization (UNCHANGED). See docs/desk/edge_mellin_symbol.md 5.2.
mpmath, high precision.
"""
import mpmath as mp
mp.mp.dps = 40

def R(m, lam):
    return (m + lam) * (m + 2) * mp.tan(mp.pi * m / 2)

def edge_blowup(lam, N=160):
    M = -lam
    a = [mp.mpf(1)]
    m = M
    for n in range(1, N + 1):
        a.append(a[-1] / R(m - 1, lam)); m -= 1
    def G(xi):
        return sum(an * xi ** (M - n) for n, an in enumerate(a))
    xis = [mp.mpf(10) ** (-k) for k in (3, 4, 5, 6)]
    xs = [1 / mp.sqrt(x) for x in xis]
    ys = [mp.log(abs(G(x))) for x in xis]
    return (ys[-1] - ys[0]) / (xs[-1] - xs[0])         # edge slope d log|G| / d(1/sqrt xi)

def winding(lam, T=2000.0, n=20000):
    taus = [(-T + 2 * T * k / n) for k in range(n + 1)]
    prev = None; tot = mp.mpf(0)
    for tau in taus:
        val = mp.arg((1j * tau + lam) * (1j * tau + 2) * (1j * mp.tanh(mp.pi * tau / 2)))
        if prev is not None:
            d = val - prev
            while d > mp.pi: d -= 2 * mp.pi
            while d < -mp.pi: d += 2 * mp.pi
            tot += d
        prev = val
    return tot / (2 * mp.pi)

print("(D) edge blow-up slope d log|G|/d(1/sqrt xi) of the terminating-up solution (lambda-indep):")
for lam in [mp.mpf('0.7'), mp.mpf('2.5'), mp.mpf('0.5')+1.5j, -mp.mpf('0.7')]:
    print(f"    lambda={str(lam):<14s} slope={mp.nstr(edge_blowup(lam),4)}  (>0 => edge-divergent)")
print("\n(W) winding W(lambda) along Re s=0  [EXACT per-factor: sign(Re lam)/2 + 1/2(s+2) + 1/2(tan)]:")
print("    NOTE: the naive line-integral is polluted by +-1/2 unwrap turns from the R(0)=0 on-line")
print("    zero; the EXACT decomposition is authoritative (+3/2 for Re lam>0, +1/2 for Re lam<0).")
for lam in [mp.mpf('0.7'), mp.mpf('2.5'), mp.mpf('-0.5')]:
    exact = (mp.sign(lam.real if hasattr(lam, 'real') else lam)) / 2 + mp.mpf(1)
    print(f"    lambda={str(lam):<8s} naive={mp.nstr(winding(lam),4)} (polluted)  EXACT={mp.nstr(exact,3)}"
          f"   (NOT the Fredholm index; tautology, see header)")
print("\n=> eigenfunctions are Mellin-line; winding is the indicial-root tautology; the leap is tau(s).")
