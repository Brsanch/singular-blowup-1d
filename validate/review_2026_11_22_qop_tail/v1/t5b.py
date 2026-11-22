import sys, numpy as np
sys.path.insert(0, '/Volumes/4TB SD/ClaudeCode/singular-blowup-1d/src')
from stability_edge_z import build_blocks_z
N, zmax = 1600, 7.
print("eigs with Re>-0.6 and |Im|<8 (Im>=0 shown), N=1600 zmax=7")
for zmin in (-10., -11., -11.25, -11.5, -12., -13., -14., -16., -18.):
    for bc, tail in ((False, True), (False, False), (True, True)):
        v = np.linalg.eigvals(build_blocks_z(N=N, zmin=zmin, zmax=zmax, edge_bc=bc, edge_tail=tail)['M_GG'])
        w = v[(v.real > -0.6) & (np.abs(v.imag) < 8) & (v.imag >= 0)]
        w = w[np.argsort(-w.real)]
        print(f" zmin={zmin:6.2f} bc={bc!s:5} tail={tail!s:5}: " + "  ".join(f"{x.real:+.4f}{x.imag:+.3f}j" for x in w[:8]) + (" ..." if len(w) > 8 else "") + ("(none)" if len(w)==0 else ""))
