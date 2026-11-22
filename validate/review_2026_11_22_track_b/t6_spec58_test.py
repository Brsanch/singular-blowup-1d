import numpy as np, math
from lib import *
from solver import *
from t6b_norms import spec58, spec14
g=Grid(-64,64,8192); u=g.u
# test: h = 1/(1+u^2):  hhat(m) = pi/sin(pi m)  on 0<Re m<1 ; add u^3 e^-u
h=1/(1+u*u)
A=spec58(g,h); xi=g.xi[:g.n//2]
ex=np.abs(Nsym(0.625+1j*xi))*np.abs(np.pi/np.sin(np.pi*(0.625+1j*np.minimum(xi,200))))
for X in (0,1,2,4,6,8,10,12):
    k=int(round(X/g.dxi)); print(X, A[k], ex[k], A[k]/ex[k]-1)
print("floor at xi=30,60:",A[int(30/g.dxi)],A[int(60/g.dxi)])
