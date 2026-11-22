import numpy as np
_x,_w=np.polynomial.legendre.leggauss(20)
def _panels(z,lo,hi):
    ed=[z]
    # left of z
    L=[];p=z;w=1e-9
    while p>lo:
        q=max(lo,p-w); L.append((q,p)); p=q; w=min(w*3,0.2)
    R=[];p=z;w=1e-9
    while p<hi:
        q=min(hi,p+w); R.append((p,q)); p=q; w=min(w*3,0.2)
    return np.array(L+R)
def J_fast(G,z,lo,hi):
    P=_panels(z,lo,hi)
    if len(P)==0: return 0.0
    a=P[:,0:1];b=P[:,1:2]
    zp=0.5*(a+b)+0.5*(b-a)*_x[None,:]
    w=0.5*(b-a)*_w[None,:]
    d=z-zp
    with np.errstate(all='ignore'):
        L=np.where(d>0, d+np.log(-np.expm1(-np.abs(d))), np.log(-np.expm1(-np.abs(d))))
    L=np.where(np.isfinite(L),L,0.0)
    return float((np.exp(zp/2)*G(zp)*L*w).sum()/np.pi)
def Q_fast(G,zs,lo,hi):
    return np.array([np.exp(-z)*J_fast(G,z,lo,hi) for z in zs])
