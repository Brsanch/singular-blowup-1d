import sys, numpy as np
from scipy.linalg import expm, svd
sys.path.insert(0, "/Volumes/4TB SD/ClaudeCode/singular-blowup-1d/src")
from stability_edge_z import build_blocks_z
for (N,zmin,zmax) in [(1400,-14.,12.),(2000,-14.,12.),(1600,-18.,12.)]:
    B=build_blocks_z(N=N,zmin=zmin,zmax=zmax,edge_bc=True); z,M=B["z"],B["M_GG"]; dz=z[1]-z[0]
    E={tau:expm(tau*M) for tau in (0.02,0.25,1.0)}
    for gam in (0.3,0.4,0.45):
        w=np.exp(-gam*z)
        # (i) where does the full-matrix top singular vector live?
        A=(w[:,None]*E[1.0])/w[None,:]
        U,S,Vh=svd(A); v=Vh[0]/w; v/=np.abs(v).max(); u=U[:,0]
        zc_in=np.sum(z*Vh[0]**2); zc_out=np.sum(z*u**2)
        rough=np.linalg.norm(np.diff(Vh[0]))/np.linalg.norm(Vh[0])
        # (ii) smooth subspace: Gaussians (width wd) centred zc in [zmin+4, 5], modulated
        for wd,lo in ((0.5,zmin+4),(0.5,zmin+1.5),(0.25,zmin+4)):
            cols=[]
            for zc in np.arange(lo,5.01,wd):
                for om in (0,1.5,3.0):
                    env=np.exp(-((z-zc)/wd)**2)
                    cols.append(env*np.cos(om*(z-zc)))
                    if om>0: cols.append(env*np.sin(om*(z-zc)))
            P=np.array(cols).T
            P0=w[:,None]*P
            Uq,Sq,Vq=svd(P0,full_matrices=False); keep=Sq>1e-9*Sq[0]; Bc=Vq[keep].T/Sq[keep]
            res=[]
            for tau,Et in E.items():
                At=(w[:,None]*(Et@P))@Bc
                # exclude the outflow strip z>zmax-3 from the output norm? keep full norm
                res.append(f"{tau}: {np.linalg.norm(At,2)*np.exp(2*gam*tau):.4f}")
            print(f"N={N} zmin={zmin} gam={gam} smooth(wd={wd},lo={lo}) dim={int(keep.sum())}: "+"  ".join(res))
        print(f"   full-matrix sv1(tau=1)*e^(2 gam)={S[0]*np.exp(2*gam):.3f}; input centroid z={zc_in:.2f}, output centroid z={zc_out:.2f}, input roughness ||diff||/||v||={rough:.3f} (2 = sawtooth)"); sys.stdout.flush()
