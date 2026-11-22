import sys, numpy as np
from scipy.linalg import expm
sys.path.insert(0, '/Volumes/4TB SD/ClaudeCode/singular-blowup-1d/src')
from stability_edge_z import build_blocks_z
TAUS = np.linspace(0.0, 6.0, 31)
for (N, zmin, zmax) in ((1200, -14, 8), (1400, -12, 8)):
    for tail in (True, False):
        B = build_blocks_z(N=N, zmin=zmin, zmax=zmax, beta=0.0, edge_bc=True, edge_tail=tail)
        M = -2.0*B["Dz"] + 0.5*B["Q_op"]
        n = [float(np.linalg.norm(expm(M*t), 2)) for t in TAUS]
        print(f" N={N} zmin={zmin} tail={tail!s:5}: worst-case sup||e^(Mt)||_2={max(n):.2f} at tau={TAUS[int(np.argmax(n))]:.1f}", flush=True)
