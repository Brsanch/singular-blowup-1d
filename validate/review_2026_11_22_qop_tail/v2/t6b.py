import sys, numpy as np
R="/Volumes/4TB SD/ClaudeCode/singular-blowup-1d"
sys.path.insert(0,R+"/validate")
import edge_coupled_h_nonlinear_closure as E
import coevolve_moving_frame as C, stability_edge_z as S
print(C.__file__, S.__file__, C.build.__defaults__)
B=E.build(N=800,zmin=-10.,zmax=10.,edge_tail=False)
B2=S.build_blocks_z(N=800,zmin=-10.,zmax=10.,edge_bc=True,edge_tail=False)
print("Qop equal to tail-off:",np.array_equal(B["Q_op"],B2["Q_op"]), "Dz eq", np.array_equal(B["Dz"],B2["Dz"]))
import inspect; print(inspect.getsource(E.check_ch4)[:900])
