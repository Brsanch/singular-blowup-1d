import sys, numpy as np
R = '/Volumes/4TB SD/ClaudeCode/singular-blowup-1d'
sys.path.insert(0, R + '/src'); sys.path.insert(0, R + '/validate')
from stability_edge_z import build_blocks_z, spectrum
def R4(x): return 1 + x + x**2/2 + x**3/6 + x**4/24
print("fig_spectrum (N=1000,zmin=-14,zmax=8,BC), window Re in (-4.5,0.4), |Im|<25")
ev = {}
for tail in (False, True):
    v, _ = spectrum(build_blocks_z(N=1000, zmin=-14.0, zmax=8.0, edge_bc=True, edge_tail=tail)["M_GG"]); ev[tail] = v
    win = (v.real > -4.5) & (v.real < 0.4) & (np.abs(v.imag) < 25)
    print(f"  tail={tail}: KEY maxRe_GG={v.real.max():.8g} npos={int((v.real>1e-6).sum())}; points in window={int(win.sum())}; min Re={v.real.min():.1f}")
a = ev[False]; b = ev[True]
wa = a[(a.real > -4.5) & (a.real < 0.4) & (np.abs(a.imag) < 25)]
d = np.array([np.abs(b - x).min() for x in wa])
print(f"  windowed truncated eigenvalues -> nearest tail eigenvalue: median shift {np.median(d):.3f}, max {d.max():.3f} (axis span 4.9 x 50)")
print("fig_gap (N=1000,zmin=-14,zmax=10,BC): rates by eigen-evolve")
taus = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
for tail in (False, True):
    B = build_blocks_z(N=1000, zmin=-14.0, zmax=10.0, edge_bc=True, edge_tail=tail); M, z = B["M_GG"], B["z"]; dz = z[1]-z[0]
    vals, V = np.linalg.eig(M); G0 = np.exp(-((z-1.0)/1.0)**2); c = np.linalg.solve(V, G0)
    Gt = [V @ (np.exp(vals*t)*c) for t in taus]
    rates = []
    for g in (0.0, 0.25, 0.5, 0.75):
        wn = lambda G: float(np.sqrt(np.sum(np.exp(-2*g*z)*np.abs(G)**2)*dz))
        r = [wn(G)/wn(G0) for G in Gt]; w = (taus >= 0.5) & (taus <= 2.5)
        rates.append(-float(np.polyfit(taus[w], np.log(np.array(r)[w]), 1)[0]))
    print(f"  tail={tail}: rates(g=0,.25,.5,.75)={np.round(rates,5)}  cond(V)={np.linalg.cond(V):.2e}")
print("uncapped RK4 (dt=0.25 dz) users: dt*|lambda_min| vs 2.785")
for nm, (N, zmin, zmax) in (("fig_coevolve build(500,-12,8)", (500,-12.,8.)), ("coevolve_projected/stress/paper_match build(1400,-12,12)", (1400,-12.,12.)), ("coevolve_projected conv (1800,-14,12)", (1800,-14.,12.)), ("coevolve_projected conv (1000,-10,12)", (1000,-10.,12.))):
    B = build_blocks_z(N=N, zmin=zmin, zmax=zmax, edge_bc=True); v = np.linalg.eigvals(B["M_GG"]); dt = 0.25*B["dz"]
    print(f"  {nm}: lambda_min={v.real.min():.1f} dt*|lambda|={dt*abs(v.real.min()):.3f} max|R4(dt lam)|={np.abs(R4(dt*v)).max():.5f}")
print("hypothetical neighbours: which (N,zmin) break uncapped dt=0.25dz, zmax=12")
for (N, zmin) in ((1800,-15.),(1400,-14.),(2400,-14.),(1800,-16.),(1000,-12.),(1000,-13.)):
    B = build_blocks_z(N=N, zmin=zmin, zmax=12., edge_bc=True); v = np.linalg.eigvals(B["M_GG"]); dt = 0.25*B["dz"]
    print(f"  N={N} zmin={zmin}: dt*|lambda|={dt*abs(v.real.min()):.3f} -> {'UNSTABLE' if np.abs(R4(dt*v)).max()>1+1e-9 else 'stable'}")
