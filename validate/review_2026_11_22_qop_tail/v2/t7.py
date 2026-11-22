import numpy as np
from common import *
OLD=load_old()
def lam(F,z,g):
    w=np.exp(-g*z); A=(w[:,None]*F)/w[None,:]; return np.linalg.eigvalsh(0.5*(A+A.T))[-1]
# beta != 0
for beta in (0.0,-0.25,-0.5,0.3):
    Bn=S.build_blocks_z(N=1000,beta=beta); B0=S.build_blocks_z(N=1000,beta=beta,edge_tail=False); Bo=OLD.build_blocks_z(N=1000,beta=beta)
    ident=all(np.array_equal(B0[k],Bo[k]) for k in Bo if isinstance(Bo[k],np.ndarray))
    v=np.linalg.eigvals(Bn['M_GG']); v0=np.linalg.eigvals(B0['M_GG'])
    # consistency: M_GG(beta) == W (-2Dz + diag(.5/ym1) J_op) W^-1 with tail
    z=Bn['z']; W=np.exp(beta*z); Mref=(W[:,None]*(-2*Bn['Dz']+np.diag(0.5*np.exp(-z))@Bn['J_op']))/W[None,:]
    print(f"beta={beta:+.2f}: tail-off==old bitwise {ident}; tail-on maxRe={v.real.max():+.5f} n+={int((v.real>1e-6).sum())} (off: {v0.real.max():+.5f}, {int((v0.real>1e-6).sum())}); |M_GG - W(T+F)W^-1|max rel={np.abs(Bn['M_GG']-Mref).max()/np.abs(Mref).max():.1e}")
# matrix-level forcing vs Q/2
for zmin in (-10.,-14.,-18.):
    B=S.build_blocks_z(N=1500,zmin=zmin); F=B['M_GG']+2*B['Dz']; d=np.abs(F-0.5*B['Q_op'])
    print(f"zmin={zmin}: max|forcing - Q_op/2| entrywise={d.max():.2e} (rel to max|Q/2| {d.max()/np.abs(0.5*B['Q_op']).max():.1e})")
# depth growth of C_force at fixed N for gamma<=0.25
for zmin in (-12.,-16.,-20.,-24.):
    B=S.build_blocks_z(N=1400,zmin=zmin,zmax=12.); 
    print(f"N=1400 zmin={zmin}: C_force(0.1)={lam(.5*B['Q_op'],B['z'],0.1):+.3e} C_force(0.2)={lam(.5*B['Q_op'],B['z'],0.2):+.3e} C_force(0.25)={lam(.5*B['Q_op'],B['z'],0.25):+.3e}")
