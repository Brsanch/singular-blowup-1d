import sys, numpy as np, pickle
d=sys.argv[1]; sys.path.insert(0,d+"/validate")
import edge_coupled_h_nonlinear_closure as E
B=E.build(N=800,zmin=-10.,zmax=10.,edge_tail=False) if 'mix2' in d else E.build(N=800,zmin=-10.,zmax=10.)
pickle.dump({k:v for k,v in B.items()},open(d+"/B.pkl","wb"))
print({k:(type(v).__name__, getattr(v,'shape',v), getattr(v,'flags',None) and (v.flags['C_CONTIGUOUS'],v.flags['F_CONTIGUOUS'])) for k,v in B.items()})
