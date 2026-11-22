import sys, numpy as np, importlib.util
def load(path,name,srcdir):
    sys.path.insert(0,srcdir)
    for m in ('stability_edge_z',): sys.modules.pop(m,None)
    sp=importlib.util.spec_from_file_location(name,path); m=importlib.util.module_from_spec(sp); sp.loader.exec_module(m)
    sys.path.pop(0); return m
R="/Volumes/4TB SD/ClaudeCode/singular-blowup-1d"; O="/private/tmp/qop_review_v2/old"
new=load(R+"/validate/coevolve_moving_frame.py","cmf_new",R+"/src")
old=load(O+"/validate/coevolve_moving_frame.py","cmf_old",O+"/src")
print(new.build_blocks_z.__module__, sys.modules.get('stability_edge_z'))
Bn=new.build(N=800,zmin=-10.,zmax=10.,edge_tail=False); Bo=old.build(N=800,zmin=-10.,zmax=10.)
for k in Bo:
    if isinstance(Bo[k],np.ndarray): print(k, np.array_equal(Bo[k],Bn[k]), np.abs(Bo[k]-Bn[k]).max())
