import numpy as np, sys, time
from common import *
from exact import *
from scipy.integrate import quad
zmax=7.0
def funcs(zmin):
    F={}
    F['bulk w/(2+w)^2']=(lambda zp: np.exp(zp)/(2+np.exp(zp))**2, -80.0, 80.0, (), True)
    F['bump c=zmin+1 r=.8']=(bump(zmin+1,.8), zmin+.2, zmin+1.8, (), False)
    F['bump c=zmin+.5 r=.4']=(bump(zmin+.5,.4), zmin+.1, zmin+.9, (), False)
    F['jump e^{zmin-z}']=(lambda zp: np.where(np.asarray(zp)>=zmin,np.exp(zmin-np.asarray(zp)),0.0)+0.0, zmin, 80.0, (), False)
    b1=bump(zmin+1,.8); b2=bump(zmin+2.5,.8)
    i1=quad(lambda t:np.exp(-t/2)*b1(t),zmin+.2,zmin+1.8)[0]; i2=quad(lambda t:np.exp(-t/2)*b2(t),zmin+1.7,zmin+3.3)[0]
    al=i1/i2
    F['zero-edge-functional b1-al*b2']=(lambda zp: b1(zp)-al*b2(zp), zmin+.2, zmin+3.3, (zmin+1.7,zmin+1.8), False)
    return F
for zmin in (-6.,-10.,-14.,-18.):
    F=funcs(zmin)
    cache={}
    for N in (400,800,1600,3200):
        Bo=NEW.build_blocks_z(N=N,zmin=zmin,zmax=zmax,edge_bc=False,edge_tail=True,hilbert_grid_support=False)
        Bn=NEW.build_blocks_z(N=N,zmin=zmin,zmax=zmax,edge_bc=False,edge_tail=True,hilbert_grid_support=True)
        z=Bo['z']
        idx=sorted(set([0,1,2,4,8]+[int(np.argmin(abs(z-t))) for t in list(zmin+np.array([.25,.5,.75,1,1.5,2,2.5,3]))+[zmin+4,zmin+5,-2.,0.,2.,4.] if t>zmin]))
        for name,(G,a,b,brk,isbulk) in F.items():
            Gv=np.array([float(G(t)) for t in z])
            ex_Q=np.array([exact_Q(G,z[i],a,b,brk) for i in idx])
            ex_H=np.array([exact_Hg(G,z[i],a,b,brk) if i>0 else np.nan for i in idx])
            if isbulk:   # also the "grid-supported" truth
                exg_Q=np.array([exact_Q(G,z[i],zmin,zmax,brk) for i in idx])
                exg_H=np.array([exact_Hg(G,z[i],zmin,zmax,brk) if (0<i<N-1) else np.nan for i in idx])
            zi=z[idx]; edge=(zi<=zmin+3.01); blk=~edge
            def err(d,e,m):
                m=m&np.isfinite(e)
                return np.nanmax(abs(d[m]-e[m]))/np.nanmax(abs(e[m])) if m.any() else np.nan
            row=[]
            for tag,B in (('old',Bo),('new',Bn)):
                dH=(B['Hg_op']@Gv)[idx]; dQ=(B['Q_op']@Gv)[idx]
                row+= [err(dH,ex_H,edge),err(dH,ex_H,blk),err(dQ,ex_Q,edge),err(dQ,ex_Q,blk)]
                if isbulk: row+=[err(dH,exg_H,edge),err(dQ,exg_Q,edge),err(dQ,exg_Q,blk)]
            print(f"zmin={zmin:5.0f} N={N:5d} {name:30s} "+"  ".join(f"{x:9.2e}" for x in row),flush=True)
print("columns: old[Hg edge, Hg bulk, Q edge, Q bulk (, vs grid-supported truth: Hg edge, Q edge, Q bulk)]  then new[same]")
print("edge = sample nodes with z<=zmin+3 (Hg excludes node 0), bulk = z>zmin+3; error = max|disc-exact|/max|exact| over the sample")
