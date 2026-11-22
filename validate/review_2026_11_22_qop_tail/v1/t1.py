import sys, importlib.util, numpy as np
def load(p, name):
    s = importlib.util.spec_from_file_location(name, p); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
old = load('/private/tmp/qop_review/old/src/stability_edge_z.py', 'old')
new = load('/Volumes/4TB SD/ClaudeCode/singular-blowup-1d/src/stability_edge_z.py', 'new')
# bit identity
allok = True
for N in (300, 801):
  for zmin in (-10., -14., -18.):
    for beta in (0.0, -0.25):
      for bc in (True, False):
        A = old.build_blocks_z(N=N, zmin=zmin, zmax=9., beta=beta, edge_bc=bc)
        B = new.build_blocks_z(N=N, zmin=zmin, zmax=9., beta=beta, edge_bc=bc, edge_tail=False)
        for k in A:
            same = np.array_equal(np.asarray(A[k]), np.asarray(B[k]))
            if not same: allok = False; print('DIFF', N, zmin, beta, bc, k)
        assert set(A) == set(B)
print('edge_tail=False bit-identical to pre-change:', allok)

# exact test functions
def fam_w(b):   # G = w/(b+w)^2 : H=(w-b)/(2 sqrt b (w+b)^2), J=(1/(2 sqrt b))[ln((w+b)/b)+2b/(w+b)-2]
    G = lambda w: w/(b+w)**2
    H = lambda w: (w-b)/(2*np.sqrt(b)*(w+b)**2)
    def Q(w):
        x = w/b
        # series for small x to avoid cancellation: ln(1+x)+2/(1+x)-2 = -x + (3/2)x^2 ... use log1p
        J = (np.log1p(x) + 2/(1+x) - 2)/(2*np.sqrt(b))
        Js = (-x + 1.5*x**2 - (5/3)*x**3 + (7/4)*x**4)/(2*np.sqrt(b))
        return np.where(x < 1e-3, Js, J)/w
    return G, H, Q
def fam_1(b):   # G = 1/(b+w)^2
    G = lambda w: 1/(b+w)**2
    H = lambda w: 1/(np.sqrt(b)*(w+b)**2) + 1/(2*b**1.5*(w+b))
    def Q(w):
        x = w/b
        J = (1/np.sqrt(b))*(1/b - 1/(w+b)) + np.log1p(x)/(2*b**1.5)
        return J/w
    return G, H, Q
G2,H2,Q2 = fam_w(2.0); G1,H1,Q1 = fam_w(1.0); c = 2**-1.5
tests = {
 'A: w/(2+w)^2 (a=-0.17678)': (G2,H2,Q2),
 'B: w/(2+w)^2 - 2^-1.5 w/(1+w)^2 (a=0, O(w))': (lambda w: G2(w)-c*G1(w), lambda w: H2(w)-c*H1(w), lambda w: Q2(w)-c*Q1(w)),
}
# check series/closed form sanity with quad
from scipy.integrate import quad
for nm,(G,H,Q) in tests.items():
    for w0 in (1e-5, 0.3, 20.):
        Jn = quad(H, 0, w0, epsabs=1e-15, epsrel=1e-13)[0]
        print('  sanity', nm[:1], w0, Q(np.array([w0]))[0], Jn/w0)
    # H closed-form vs numeric pv at w0=0.7
    w0=0.7
    f = lambda s: 2*s*s*G(s*s)/(s*s)  # g(t)=t^-1/2 G(t), t=s^2, dt=2s ds -> 2 G(s^2)
    from scipy.integrate import quad as q
    s0=np.sqrt(w0)
    # pv int_0^inf 2G(s^2)/(w0-s^2) ds = pv int 2G/((s0-s)(s0+s))
    h = lambda s: -2*G(s*s)/(s0+s)
    I1 = q(h, 0, 2*s0, weight='cauchy', wvar=s0)[0]
    I2 = q(lambda s: 2*G(s*s)/(w0-s*s), 2*s0, np.inf)[0]
    print('  H sanity', nm[:1], H(np.array([w0]))[0], (I1+I2)/np.pi)

zmax = 12.0
for nm,(G,H,Q) in tests.items():
    print('\n===', nm)
    for zmin in (-10., -14., -18.):
        for N in (1500, 3000):
            row = []
            for bc in (True, False):
              for tail in (False, True):
                B = new.build_blocks_z(N=N, zmin=zmin, zmax=zmax, edge_bc=bc, edge_tail=tail)
                z = B['z']; w = np.exp(z)
                Gv = G(w); Qe = Q(w); Qn = B['Q_op'] @ Gv
                He = H(w); Hn = B['Hg_op'] @ Gv
                # forcing lives inside M_GG = -2Dz + forcing ; recover
                F = B['M_GG'] + 2*B['Dz']
                Fn = 2*(F @ Gv)
                err = np.abs(Qn-Qe); m = z < zmax-4   # exclude far-field truncation zone separately
                i = np.argmax(err[m])
                rel_edge = (Qn[0]-Qe[0])/Qe[0]
                row.append((bc, tail, err[m].max(), z[m][i], Qn[0], Qe[0], rel_edge, np.abs(Fn-Qn)[m].max(), np.abs(Fn-Qe)[m].max(), np.abs(Hn-He)[m].max(), err[~m].max()))
            for r in row:
                print(f" zmin={zmin:5.0f} N={N} bc={r[0]!s:5} tail={r[1]!s:5} max|Qop G-Q|={r[2]:.3e} at z={r[3]:+.2f}  Qn[first]={r[4]:+.4e} Qex={r[5]:+.4e} relerr_first={r[6]:+.3f}  max|2F G-Qop G|={r[7]:.2e} max|2F G-Q|={r[8]:.2e} max|Hg-H|={r[9]:.2e} farzone={r[10]:.1e}")
