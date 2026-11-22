import sys, numpy as np
sys.path.insert(0, "/Volumes/4TB SD/ClaudeCode/singular-blowup-1d/src")
import stability_edge_z as S
b = 2.0
for (N, zmin, zmax) in ((1400, -14, 7), (2800, -14, 7), (2800, -20, 12), (4000, -24, 14)):
    B = S.build_blocks_z(N=N, zmin=zmin, zmax=zmax, edge_bc=False)
    z = B["z"]; w = np.exp(z)
    G = w / (b + w) ** 2
    qex = np.where(w < 1e-8, -1 / (2 * b**1.5) + 3 * w / (4 * b**2.5), np.log1p(w / b) / (2 * np.sqrt(b) * w) - 1 / (np.sqrt(b) * (b + w)))
    q = B["Q_op"] @ G; Hh = B["Hg_op"] @ G
    m = (z > zmin + 4) & (z < zmax - 4)
    print(f"repo Q_op N={N} z=[{zmin},{zmax}]: max|Q_op G - Q_exact| interior={np.max(np.abs(q-qex)[m]):.2e}; Hg at edge={Hh[0]:+.6f}, Q_op row1={q[1]:+.6f} (exact a={-1/(2*b**1.5):+.6f}); Q_op[0]={q[0]:+.3e}")
