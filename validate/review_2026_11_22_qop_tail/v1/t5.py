import sys, numpy as np
sys.path.insert(0, '/Volumes/4TB SD/ClaudeCode/singular-blowup-1d/src')
from stability_edge_z import build_blocks_z
np.set_printoptions(linewidth=250)
N, zmax = 1600, 7.
print("free edge + tail: all eigs with Re>-0.45, |Im|<8  (N=1600,zmax=7); and BC+tail same window")
for zmin in (-10., -11., -11.25, -11.5, -11.75, -12., -12.5, -13., -14., -16., -18.):
    B = build_blocks_z(N=N, zmin=zmin, zmax=zmax, edge_bc=False, edge_tail=True)
    v, V = np.linalg.eig(B['M_GG']); z = B['z']
    m = (v.real > -0.45) & (np.abs(v.imag) < 8) & (v.imag >= 0)
    sel = np.argsort(-v.real * m - 1e9*(~m))[:int(m.sum())]
    s = []
    for i in sel:
        vec = np.abs(V[:, i]); k = (z < zmin + 4) & (z > zmin + 0.3)
        alpha = np.polyfit(z[k], np.log(vec[k] + 1e-300), 1)[0]
        s.append(f"{v[i].real:+.4f}{v[i].imag:+.3f}j(edge slope {alpha:+.2f}; -Re/2={-v[i].real/2:+.2f})")
    vb = np.linalg.eigvals(build_blocks_z(N=N, zmin=zmin, zmax=zmax, edge_bc=True, edge_tail=True)['M_GG'])
    mb = (vb.real > -0.45) & (np.abs(vb.imag) < 8)
    print(f" zmin={zmin:6.2f} free+tail: n={int(m.sum())} " + "  ".join(s) + f" || BC+tail in window: n={int(mb.sum())}")
print("\nN-dependence of the crossing: max Re of free+tail")
for N in (600, 1000, 1600, 2400):
    r = []
    for zmin in (-11., -11.25, -11.5, -11.75, -12.):
        v = np.linalg.eigvals(build_blocks_z(N=N, zmin=zmin, zmax=zmax, edge_bc=False, edge_tail=True)['M_GG'])
        r.append(f"{zmin}:{v.real.max():+.4f}")
    print(f" N={N}: " + "  ".join(r))
