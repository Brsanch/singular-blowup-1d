import sys, numpy as np
sys.path.insert(0,"/private/tmp/qop_review_v2/mix2/validate")
import edge_coupled_h_nonlinear_closure as E, stability_edge_z as S
orig=S.build_blocks_z
calls=[]
def spy(*a,**k):
    calls.append(k); return orig(*a,**k)
import coevolve_moving_frame as C
C.build_blocks_z=spy
E.check_ch4(); print(calls)
