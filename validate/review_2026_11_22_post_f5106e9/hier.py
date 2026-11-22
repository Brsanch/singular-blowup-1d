# Review script (2026-11-22, post-f5106e9 review): min eigenvalue of the pointwise hierarchical u-dissipation form, k = 2, 3, 4.
import numpy as np, math
def lam(c):
    k=len(c)-1; A=np.zeros((k+1,k+1))
    for m in range(k+1):
        for j in range(m+1): A[m,j]=c[m]*math.comb(m,j)*(-0.5)**(m-j)
    S=0.5*(A+A.T); D=np.diag(1/np.sqrt(c)); return np.linalg.eigvalsh(D@S@D).min()
for k,base in [(2,4),(3,4),(4,6),(4,4),(3,6),(4,1),(3,1)]:
    c=np.array([float(base)**(k-m) for m in range(k+1)]); print(k,base,round(lam(c),6))
# IBP-improved bound: integrate Re v_j conj v_m parts? just report pointwise
