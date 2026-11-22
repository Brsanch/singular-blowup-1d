import numpy as np, sys
from multiprocessing import Pool
from common import *; from fastq import *
OLD=load_old()
def run(zmin):
    out=[]
    zmax=7.
    fs={'a':Ga,'b':mk_bump(zmin+1),'c':mk_jump(zmin),'d':mk_azero(zmin),'e':mk_bump(1.0,1.0)}
    for N in (800,1600,3200):
        zfull=np.linspace(zmin,zmax,N)
        ref={}
        for k,G in fs.items():
            ref[k,'grid']=Q_fast(G,zfull,zmin,zmax)
            ref[k,'cont']=Q_fast(G,zfull,-70.,45.) if k!='c' else None
        for bc in (False,True):
            B=S.build_blocks_z(N=N,zmin=zmin,zmax=zmax,edge_bc=bc)
            B0=S.build_blocks_z(N=N,zmin=zmin,zmax=zmax,edge_bc=bc,edge_tail=False)
            Bo=OLD.build_blocks_z(N=N,zmin=zmin,zmax=zmax,edge_bc=bc)
            ident=all(np.array_equal(B0[k],Bo[k]) for k in Bo if isinstance(Bo[k],np.ndarray))
            z=B['z']; sl=slice(1,None) if bc else slice(None)
            for k,G in fs.items():
                Gv=G(z)
                q=B['Q_op']@Gv; q0=B0['Q_op']@Gv
                f2=np.exp(-z)*(B['J_op']@Gv)
                Fm=B['M_GG']+2*B['Dz']
                f3=2*(Fm@Gv)
                sc=np.abs(ref[k,'grid'][sl]).max()
                eg=np.abs(q-ref[k,'grid'][sl]); 
                ec=np.abs(q-ref[k,'cont'][sl]).max() if ref[k,'cont'] is not None else np.nan
                e0=np.abs(q0-ref[k,'grid'][sl]).max()
                relpt=(eg/np.maximum(np.abs(ref[k,'grid'][sl]),1e-3*sc)).max()
                out.append(f"zmin={zmin:5.0f} N={N:4d} bc={int(bc)} G={k} oldident={ident} scale={sc:9.3e} maxerr_grid={eg.max():9.3e} (at z={z[eg.argmax()]:+.2f}) rel={eg.max()/sc:8.2e} maxptrel={relpt:8.2e} maxerr_cont={ec:9.3e} truncated_err={e0:9.3e} | Q-e^-zJ={np.abs(q-f2).max():8.1e} Q-2forcing={np.abs(q-f3).max():8.1e} edgeval={q[0]:+.6f} ref={ref[k,'grid'][sl][0]:+.6f}")
    return out
if __name__=="__main__":
    with Pool(3) as p:
        for o in p.map(run,(-10.,-14.,-18.)):
            print("\n".join(o)); sys.stdout.flush()
