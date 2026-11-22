"""Task 6 (a): validate reviewer's solver: temporal order, grid/domain convergence, and the
REAL-SPACE residual of  G_s - (1/2) Q G - P(G) = 0  (real-space Q, s-derivative from snapshots)."""
import math, sys, time
import numpy as np
from lib import *
from solver import *

A = float(sys.argv[1]) if len(sys.argv) > 1 else 0.05   # max|G_0| amplitude (gate convention)
S = 1.5
sc = G0_maxnorm_scale()
print(f"amplitude A=max|G0|={A}; max|N h0|={sc:.6f}")


def solve(zmin, zmax, n, nsteps, store_every):
    s = Solver(zmin, zmax, n, sponge=10.0)
    h0 = A * h0_family(s) / sc
    t0 = time.time()
    H, bad = s.run(h0, S, nsteps, store_every)
    return s, H, bad, time.time() - t0


# temporal convergence at fixed grid
ref = None
sols = {}
for nst in (150, 300, 600, 1200):
    s, H, bad, dt = solve(-64.0, 64.0, 4096, nst, nst // 150)
    G = np.array([s.N(h) for h in H])
    sols[nst] = (s, H, G)
    print(f" nsteps={nst}: {dt:.1f}s nan={bad}")
bulk = (sols[150][0].g.z > -12) & (sols[150][0].g.z < 12)
edge = (sols[150][0].g.z > -18) & (sols[150][0].g.z < -12)
for a_, b_ in ((150, 300), (300, 600), (600, 1200)):
    dG = np.abs(sols[a_][2] - sols[b_][2])
    w = sols[150][0].g.w
    print(f" |G_{a_} - G_{b_}|: bulk {dG[:, bulk].max():.3e}   edge(G/w) {(dG/w)[:, edge].max():.3e}")

# grid / domain convergence at nsteps=600
s1, H1, G1 = sols[600]
for (zmin, zmax, n) in ((-64.0, 64.0, 8192), (-96.0, 96.0, 8192), (-104.0, 24.0, 8192)):
    s2, H2, bad, dt = solve(zmin, zmax, n, 600, 4)
    G2 = np.array([s2.N(h) for h in H2])
    # compare on common points by spectral-quality interpolation (cubic is enough for a check)
    from scipy.interpolate import CubicSpline
    errb = 0.0
    for j in (0, 50, 100, 150):
        cs = CubicSpline(s2.g.z, G2[j])
        errb = max(errb, float(np.max(np.abs(cs(s1.g.z[bulk]) - G1[j][bulk]))))
    print(f" grid [{zmin},{zmax}] n={n}: nan={bad} ({dt:.1f}s)  max bulk |G - G_ref| = {errb:.3e}")

# real-space residual at nsteps=1200 (store every 8 -> snapshot spacing 0.01), 4th-order FD in s
s, H, G = sols[1200]
z = s.g.z
dsn = S / 150
Aq = Q_real_matrix(z)
print("\n REAL-SPACE residual  R = G_s - (1/2) Q_real G - P_real(G)   (per region max |R|; ref scales)")
print("   s      edge(-30,-12)   core(-12,0)    bulk(0,12)    far(12,30)  | max|G_s|   max|NL|   max|(1/2)QG|   max|Qreal-Qconj|   R with NL dropped (core)")
for j in (10, 40, 70, 100, 130, 148):
    Gs = (-G[j + 2] + 8 * G[j + 1] - 8 * G[j - 1] + G[j - 2]) / (12 * dsn)
    Qr = Q_real(G[j], z, Aq)
    DG = np.gradient(G[j], s.g.dz, edge_order=2)  # FD derivative: no spectral operator in the residual
    DGs = s.D(G[j])
    NLfd = -Qr * DG + 0.5 * Qr * G[j]
    NL = -Qr * DGs + 0.5 * Qr * G[j]
    R = Gs - 0.5 * Qr - NL
    Qc = s.Q_of_h(H[j])
    reg = [((z > lo) & (z < hi)) for lo, hi in ((-30, -12), (-12, 0), (0, 12), (12, 30))]
    core = reg[1] | reg[2]
    print(f"  {j*dsn:5.2f}   " + "   ".join(f"{np.max(np.abs(R[m])):.2e}" for m in reg)
          + f"  | {np.max(np.abs(Gs[core])):.2e}  {np.max(np.abs(NL[core])):.2e}  {np.max(np.abs(0.5*Qr[core])):.2e}   {np.max(np.abs(Qr-Qc)[core]):.2e}"
          + f"          {np.max(np.abs((Gs-0.5*Qr)[reg[1]])):.2e}   [FD-D variant core: {np.max(np.abs((Gs-0.5*Qr-NLfd)[core])):.2e}]")
