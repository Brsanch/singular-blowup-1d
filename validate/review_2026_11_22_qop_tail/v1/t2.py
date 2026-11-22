import sys, numpy as np
sys.path.insert(0, '/Volumes/4TB SD/ClaudeCode/singular-blowup-1d/src')
from stability_edge_z import build_blocks_z
np.set_printoptions(linewidth=200, precision=5, suppress=False)
def R4(x): return 1 + x + x**2/2 + x**3/6 + x**4/24
def hand(N, zmin, zmax):
    # hand-tail operator used by the three gates: -2Dz + 0.5*(Q_op_trunc + e^{zmin-z} (x) h1)
    B = build_blocks_z(N=N, zmin=zmin, zmax=zmax, edge_bc=True, edge_tail=False)
    z, dz = B['z'], B['dz']; w = np.full(len(z), dz); w[0]*=.5; w[-1]*=.5
    h1 = -(1/np.pi)*np.exp(-0.5*z)*w
    return -2*B['Dz'] + 0.5*(B['Q_op'] + np.outer(np.exp(zmin - z), h1)), h1
print("=== (a) spectra: rightmost 6, leftmost real, rho; four combos + hand-tail(BC)")
for (N, zmin, zmax) in [(1000,-12.,7.), (1600,-12.,7.), (800,-14.,7.), (1600,-14.,7.), (2400,-14.,7.), (1600,-18.,7.), (1600,-10.,7.)]:
    print(f"\n N={N} zmin={zmin} zmax={zmax}  dz={(zmax-zmin)/(N-1):.5f}  e^(-zmin/2)={np.exp(-zmin/2):.1f}")
    ev = {}
    for bc in (False, True):
        for tail in (False, True):
            B = build_blocks_z(N=N, zmin=zmin, zmax=zmax, edge_bc=bc, edge_tail=tail)
            v, V = np.linalg.eig(B['M_GG']); ev[bc,tail] = v
            o = np.argsort(-v.real)
            il = np.argmin(v.real); 
            vec = V[:, il]; vec = vec/np.abs(vec).max()
            loc = np.abs(vec[:4])
            print(f"  bc={bc!s:5} tail={tail!s:5} top: " + " ".join(f"{x.real:+.5f}{x.imag:+.3f}j" for x in v[o[:6]]) +
                  f" | n(Re>-0.01)={int((v.real>-0.01).sum())} n(Re>0)={int((v.real>1e-9).sum())} | min Re eig={v[il].real:+.2f}{v[il].imag:+.2f}j  rho={np.abs(v).max():.2f} max|Im|={np.abs(v.imag).max():.2f} | |evec| first4 nodes={loc}  mass beyond node 10={np.abs(vec[10:]).max():.2e}")
    L, h1 = hand(N, zmin, zmax); v = np.linalg.eigvals(L); il=np.argmin(v.real)
    print(f"  hand-tail (BC, h1=H[g](1) trapezoid): max Re={v.real.max():+.5f} min Re eig={v[il].real:+.2f} rho={np.abs(v).max():.2f}   ratio to e^(-zmin/2): {v[il].real/np.exp(-zmin/2):+.4f}")
    # how different is free+tail spectrum from BC+tail spectrum (set distance on the right part)
    a = ev[False,True]; b = ev[True,True]
    ar = a[a.real > -1.0]; br = b[b.real > -1.0]
    d = np.array([np.abs(br - x).min() for x in ar])
    print(f"  free+tail vs BC+tail, eigs with Re>-1: counts {len(ar)} vs {len(br)}; max nearest-distance {d.max():.2e}")
    a = ev[True,False]
    d2 = np.array([np.abs(a - x).min() for x in br])
    print(f"  BC+tail vs BC+trunc, eigs with Re>-1: counts {len(br)} vs {int((a.real>-1).sum())}; max nearest-distance {d2.max():.2e}")
    lam = ev[True,True]; il = np.argmin(lam.real)
    print(f"  BC+tail big eig / e^(-zmin/2) = {lam[il].real/np.exp(-zmin/2):+.4f};  * dz = {lam[il].real*(zmax-zmin)/(N-1):+.3f}")
