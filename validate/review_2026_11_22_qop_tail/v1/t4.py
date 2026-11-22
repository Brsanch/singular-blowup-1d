import sys, numpy as np
sys.path.insert(0, '/Volumes/4TB SD/ClaudeCode/singular-blowup-1d/src')
from stability_edge_z import build_blocks_z
def R4(x): return 1 + x + x**2/2 + x**3/6 + x**4/24
print("=== (b) author's eigenvalue claims  (-106/-702/-1842 at zmin -10/-14/-16)")
for (N,zmin,zmax) in [(1400,-10.,12.),(1400,-14.,12.),(1400,-16.,12.),(1800,-16.,14.),(1100,-12.,12.),(700,-14.,12.),(2800,-14.,12.)]:
    B = build_blocks_z(N=N, zmin=zmin, zmax=zmax, edge_bc=True, edge_tail=True)
    v, V = np.linalg.eig(B['M_GG']); il = np.argmin(v.real); dz = B['dz']
    B0 = build_blocks_z(N=N, zmin=zmin, zmax=zmax, edge_bc=True, edge_tail=False)
    v0 = np.linalg.eigvals(B0['M_GG'])
    # rank-one prediction: 0.5 * Hg_op_full[0,1:] @ e^{zmin - z}
    Bf = build_blocks_z(N=N, zmin=zmin, zmax=zmax, edge_bc=False, edge_tail=False)
    pred = 0.5 * Bf['Hg_op'][0,1:] @ np.exp(zmin - B['z'])
    vec = (V[:, il]/V[0, il]).real; z = B['z']
    shape = np.exp(zmin - z); shape = shape/shape[0]
    fitq = np.abs(vec - shape).max()
    print(f" N={N} zmin={zmin} zmax={zmax} dz={dz:.5f}: min Re eig={v[il].real:+.2f}{v[il].imag:+.1e}j  /e^(-zmin/2)={v[il].real/np.exp(-zmin/2):+.4f}  rank-one pred 0.5*Hg[0]@e^(zmin-z)={pred:+.2f}  -(1/2pi)e^(-zmin/2)(ln(1/dz)+C): C={-v[il].real*2*np.pi/np.exp(-zmin/2)-np.log(1/dz):+.3f} | rho tail={np.abs(v).max():.1f} trunc={np.abs(v0).max():.1f} | evec vs e^(zmin-z): max dev={fitq:.3f}; evec at z-zmin=0.5,1,2,4: {[round(float(np.interp(zmin+d, z, vec)),4) for d in (0.5,1,2,4)]} (e^-d: {[round(float(np.exp(-d)),4) for d in (0.5,1,2,4)]})")
    for name, dt in (("uncapped .25dz", 0.25*dz), ("capped min(.25dz,2/rho)", min(0.25*dz, 2/np.abs(v).max()))):
        g = np.abs(R4(dt*v)); gh = np.abs(R4(dt*np.linalg.eigvals(B['M_hh']))) if name.startswith("capped") else np.array([0.])
        print(f"     {name}: dt={dt:.5f} dt*rho={dt*np.abs(v).max():.3f} max|R4(dt*lam)| over eig(M_GG)={g.max():.6f} (n>1+1e-12: {int((g>1+1e-12).sum())})  M_hh: {gh.max():.6f}")
print("\n=== (d) coevolve_moving_frame grids (RK4, dt=0.25 dz, NOT capped)")
for (N,zmin,zmax) in [(1000,-10.,12.),(1400,-12.,12.),(1800,-14.,12.)]:
    B = build_blocks_z(N=N, zmin=zmin, zmax=zmax, edge_bc=True, edge_tail=True); v = np.linalg.eigvals(B['M_GG']); dt=0.25*B['dz']
    print(f" N={N} zmin={zmin}: min Re eig={v.real.min():+.1f} dt*|lam|={dt*np.abs(v.real.min()):.3f} (limit 2.785) max|R4|={np.abs(R4(dt*v)).max():.5f}")
print("\n=== (e) matrix-level: forcing vs 0.5*Q_op ; (f) hand tail vs module tail")
for zmin in (-10., -14., -18.):
  for bc in (True, False):
    for tail in (False, True):
        B = build_blocks_z(N=1000, zmin=zmin, zmax=10., edge_bc=bc, edge_tail=tail)
        F = B['M_GG'] + 2*B['Dz']; D = F - 0.5*B['Q_op']
        print(f" zmin={zmin} bc={bc!s:5} tail={tail!s:5}: max|forcing-0.5Q_op|={np.abs(D).max():.2e} max|forcing|={np.abs(F).max():.2e} rel={np.abs(D).max()/np.abs(F).max():.1e}  J_op vs diag(ym1)Q_op rel={np.abs(B['J_op']-np.exp(B['z'])[:,None]*B['Q_op']).max()/np.abs(B['J_op']).max():.1e}")
  B0 = build_blocks_z(N=1000, zmin=zmin, zmax=10., edge_bc=True, edge_tail=False)
  B1 = build_blocks_z(N=1000, zmin=zmin, zmax=10., edge_bc=True, edge_tail=True)
  Bf = build_blocks_z(N=1000, zmin=zmin, zmax=10., edge_bc=False, edge_tail=False)
  z, dz = B0['z'], B0['dz']; w = np.full(len(z), dz); w[0]*=.5; w[-1]*=.5
  h1 = -(1/np.pi)*np.exp(-0.5*z)*w
  hm = Bf['Hg_op'][0,1:]
  # module tail row recovered
  rowm = (B1['Q_op'] - B0['Q_op'])[0]*np.exp(z[0]-zmin)
  print(f"  module tail row == Hg_op_full[0,1:]: {np.abs(rowm-hm).max()/np.abs(hm).max():.1e};  hand row h1 vs module row: first 3 entries h1={h1[:3]} module={hm[:3]}; ||h1-hm||_2/||h1||_2={np.linalg.norm(h1-hm)/np.linalg.norm(h1):.2f}")
  for nm, G in (("w/(2+w)^2", np.exp(z)/(2+np.exp(z))**2), ("gauss(z-1)", np.exp(-(z-1)**2)), ("gauss((z-zmin-1)/0.3) near-edge", np.exp(-((z-zmin-1)/0.3)**2)), ("e^(zmin-z) jump mode", np.exp(zmin-z))):
      print(f"    a functional on {nm:34s}: hand h1@G={h1@G:+.5e}  module Hg[0]@G={hm@G:+.5e}")
