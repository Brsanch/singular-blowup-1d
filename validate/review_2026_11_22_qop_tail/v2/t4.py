import numpy as np, sys
from multiprocessing import Pool
from common import *
def one(a):
    N,zmin,zmax=a; s=f"N={N:4d} zmin={zmin:5.1f} zmax={zmax:4.1f} |"
    dz=(zmax-zmin)/(N-1)
    for bc,tail in ((False,False),(False,True),(True,False),(True,True)):
        v=np.linalg.eigvals(S.build_blocks_z(N=N,zmin=zmin,zmax=zmax,edge_bc=bc,edge_tail=tail)['M_GG'])
        low=v[np.abs(v.imag)<20]
        s+=f" bc{int(bc)}t{int(tail)}: all {v.real.max():+.5f} low {low.real.max():+.3f} n+={int((v.real>1e-9).sum())}"
        if (not bc) and tail:
            # anything to the right of the low-frequency maximum, outside the window?
            lm=low.real.max()
            out=v[(np.abs(v.imag)>=20)&(v.real>lm)]
            saw=np.abs(np.abs(out.imag)-2/dz)<0.15*2/dz
            mid=out[~saw]
            extra=f" [free+tail: {len(out)} eigs right of low-max outside window; non-sawtooth(|Im| not within 15% of 2/dz): {len(mid)}"
            if len(mid): 
                j=np.argmax(mid.real); extra+=f", rightmost {mid[j].real:+.4f}{mid[j].imag:+.1f}i"
            j=np.argmax(v.real); extra+=f"; overall rightmost {v[j].real:+.5f}{v[j].imag:+.1f}i (2/dz={2/dz:.1f})]"
            s+=extra
    return s
if __name__=="__main__":
    cfg=[(1600,float(z),7.) for z in (-6,-10,-12,-14,-16,-18)]+[(800,float(z),7.) for z in range(-18,-3)]+[(800,-11.5,7.),(800,-5.5,7.),(800,-4.5,7.)]+[(1200,float(z),7.) for z in (-4,-6,-10,-14,-18)]+[(N,float(z),7.) for N in (100,200,400) for z in (-4,-6,-10,-14)]+[(800,float(z),zm) for z in (-4,-6,-10,-14) for zm in (3.,12.,16.)]
    with Pool(4) as p:
        for o in p.imap(one,cfg): print(o); sys.stdout.flush()
