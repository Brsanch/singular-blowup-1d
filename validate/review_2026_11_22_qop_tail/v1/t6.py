import sys, numpy as np
from scipy.integrate import quad, dblquad
sys.path.insert(0, '/Volumes/4TB SD/ClaudeCode/singular-blowup-1d/src')
from stability_edge_z import build_blocks_z
# exact tail functional for G supported on z>zmin:  T[G] = (1/wmin) int_0^wmin H(g)(w) dw
#   = -(1/pi) int e^{z'/2} G(z') e^{-zmin} * (-ln(1 - e^{zmin - z'})) dz'
for zmin in (-10., -14.):
    wmin = np.exp(zmin)
    for nm, Gf in (("near-edge gauss (c=zmin+1, s=0.3)", lambda z: np.exp(-((z-zmin-1)/0.3)**2)),
                   ("bulk gauss (c=1, s=1)", lambda z: np.exp(-(z-1)**2)),
                   ("jump mode e^(zmin-z)", lambda z: np.exp(zmin-z))):
        Tk = -(1/np.pi)*quad(lambda s: np.exp((zmin+s)/2)*Gf(zmin+s)*np.exp(-zmin)*(-np.log1p(-np.exp(-s))), 0, 40, limit=400, points=[1.0])[0]
        # brute force: average over w in (0,wmin) of H(g)(w) -- only for the first (smooth, compact) case
        if nm.startswith("near"):
            Hw = lambda w: (1/np.pi)*quad(lambda zp: np.exp(zp/2)*Gf(zp)/(w-np.exp(zp)), zmin, zmin+4, limit=400, points=[zmin+1])[0]
            Tb = quad(lambda u: Hw(u*wmin), 0, 1, limit=200)[0]
        else: Tb = float('nan')
        out = []
        for N in (1000, 2000):
            B = build_blocks_z(N=N, zmin=zmin, zmax=10., edge_bc=True, edge_tail=False); Bf = build_blocks_z(N=N, zmin=zmin, zmax=10., edge_bc=False, edge_tail=False)
            z, dz = B['z'], B['dz']; w = np.full(len(z), dz); w[0]*=.5; w[-1]*=.5
            h1 = -(1/np.pi)*np.exp(-0.5*z)*w; hm = Bf['Hg_op'][0,1:]
            out.append(f"N={N}: hand={h1@Gf(z):+.4f} module={hm@Gf(z):+.4f}")
        print(f"zmin={zmin} {nm:36s}: exact kernel T={Tk:+.4f} (brute {Tb:+.4f}) | " + " ; ".join(out) + f" | T/e^(-zmin/2)={Tk/np.exp(-zmin/2):+.4f}")
    # eigenvalue of BC operator with the exact-kernel tail row (midpoint-ish quadrature of the log kernel)
    for N in (1000, 2000):
        B = build_blocks_z(N=N, zmin=zmin, zmax=10., edge_bc=True, edge_tail=False); z, dz = B['z'], B['dz']
        # integrate log kernel exactly over each cell [z_j - dz/2, z_j + dz/2] (first cell from zmin): weight_j = int k(s) ds, k(s) = -ln(1-e^{-s}) e^{s/2}
        edges = np.concatenate(([zmin], 0.5*(z[1:]+z[:-1]), [z[-1]])) - zmin
        kw = np.array([quad(lambda s: -np.log1p(-np.exp(-s))*np.exp(s/2), a, b)[0] if a < 30 else 0.0 for a, b in zip(edges[:-1], edges[1:])])
        tau = -(1/np.pi)*np.exp(-zmin/2)*kw
        L = -2*B['Dz'] + 0.5*(B['Q_op'] + np.outer(np.exp(zmin - z), tau))
        v = np.linalg.eigvals(L)
        print(f"   exact-kernel-tail BC operator N={N}: min Re eig={v.real.min():+.2f} (/e^(-zmin/2)={v.real.min()/np.exp(-zmin/2):+.4f}; rank-one continuum prediction -(1/2pi)(4-4ln2)={-(4-4*np.log(2))/(2*np.pi):+.4f}), max Re={v.real.max():+.5f}")
