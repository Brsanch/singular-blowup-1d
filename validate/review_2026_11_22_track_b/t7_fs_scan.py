"""FS-3 solver (gate code, unmodified numerics) on a (domain, Nz, eps) scan: finite? and H3 / a_quad."""
import io, contextlib, warnings, re
from mut import SRC
def one(zmin,zmax,nz,eps=0.03,sponge=False):
    s=SRC["FS"]
    reps=[("        zmin, zmax = -72.0, 72.0\n",f"        zmin, zmax = {zmin}, {zmax}\n"),
          ("for Nz in (4096, 8192):",f"for Nz in ({nz},):"),
          ("    drifts = [abs(results[8192][i] - results[4096][i]) / results[8192][i] for i in range(4)]\n    a_q = abs(results[8192][4])","    drifts = [0.0]\n    a_q = abs(results[%d][4])"%nz),
          ("results[8192]","results[%d]"%nz),
          ("    eps = 0.03\n    S_max = 1.2\n    Ns = 240\n    ds = S_max / Ns\n\n    results",f"    eps = {eps}\n    S_max = 1.2\n    Ns = 240\n    ds = S_max / Ns\n\n    results")]
    if sponge:
        reps.append(("            return NgI(-QG_cur * Dg(G_cur) + 0.5 * QG_cur * G_cur)","            return NgI(-QG_cur * Dg(G_cur) + 0.5 * QG_cur * G_cur) * (0.5*(1+np.tanh((z-z[0]-8.0)/1.5))) * (0.5*(1+np.tanh((z[-1]-8.0-z)/1.5)))"))
    for o,n in reps:
        assert o in s, o
        s=s.replace(o,n)
    ns={"__name__":"m"}; buf=io.StringIO()
    with contextlib.redirect_stdout(buf), warnings.catch_warnings():
        warnings.simplefilter("ignore"); exec(compile(s,"m","exec"),ns); ns["check_fs3"]()
    t=buf.getvalue()
    m=re.search(r"H2=(\S+), H3=(\S+) .*a_quad\|=(\S+)",t)
    return f"H2={m.group(1)} H3={m.group(2)} a={m.group(3)}"
for sp in (False,True):
    print("=== end-sponge on forcing:",sp)
    for eps in (0.03,):
        for dom in ((-36,36),(-48,48),(-64,64),(-72,72),(-80,80),(-96,96),(-104,24),(-128,128)):
            row=[]
            for nz in (2048,4096,8192,16384):
                row.append(f"Nz={nz}: {one(float(dom[0]),float(dom[1]),nz,eps,sp)}")
            print(f" eps={eps} [{dom[0]},{dom[1]}]  "+" | ".join(row),flush=True)
