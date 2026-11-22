import sys, numpy as np
sys.path.insert(0, '/Volumes/4TB SD/ClaudeCode/singular-blowup-1d/src')
from stability_edge_z import build_blocks_z
def R4(x): return 1 + x + x**2/2 + x**3/6 + x**4/24
print("=== (a) free edge: real eigenvalues & low-frequency rightmost, vs zmin")
for (N, zmax) in [(800, 7.), (1600, 7.), (1200, 12.)]:
    for zmin in (-6., -8., -9., -10., -11., -11.5, -12., -13., -14., -16., -18.):
        out = []
        for tail in (False, True):
            v = np.linalg.eigvals(build_blocks_z(N=N, zmin=zmin, zmax=zmax, edge_bc=False, edge_tail=tail)['M_GG'])
            re = np.sort(v[np.abs(v.imag) < 1e-8].real)[::-1][:3]
            lf = v[np.abs(v.imag) < 20]; lfm = lf[np.argmax(lf.real)]
            out.append(f"tail={tail!s:5} maxRe={v.real.max():+.5f} top real eigs={np.round(re,4)} rightmost |Im|<20: {lfm.real:+.4f}{lfm.imag:+.3f}j")
        vb = np.linalg.eigvals(build_blocks_z(N=N, zmin=zmin, zmax=zmax, edge_bc=True, edge_tail=True)['M_GG'])
        lf = vb[np.abs(vb.imag) < 20]; lfm = lf[np.argmax(lf.real)]
        vb0 = np.linalg.eigvals(build_blocks_z(N=N, zmin=zmin, zmax=zmax, edge_bc=True, edge_tail=False)['M_GG'])
        lf0 = vb0[np.abs(vb0.imag) < 20]; lfm0 = lf0[np.argmax(lf0.real)]
        print(f" N={N} zmax={zmax} zmin={zmin:6.1f} | " + " | ".join(out) + f" | BC+tail rightmost |Im|<20: {lfm.real:+.4f}{lfm.imag:+.3f}j ; BC+trunc: {lfm0.real:+.4f}{lfm0.imag:+.3f}j")
