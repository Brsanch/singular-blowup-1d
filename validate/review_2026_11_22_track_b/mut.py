import io, contextlib, time, warnings
V="/Volumes/4TB SD/ClaudeCode/singular-blowup-1d/validate/"
SRC={k:open(V+f).read() for k,f in (("PM","edge_pure_step_nonlinear_manifold.py"),("NM","edge_nonlinear_cusp_manifold.py"),("FS","edge_pure_step_finite_sobolev.py"))}
def run(label, key, fn, reps):
    s=SRC[key]
    for old,new in reps:
        assert old in s, f"pattern not found: {old!r}"
        s=s.replace(old,new)
    ns={"__name__":"mut"}; buf=io.StringIO(); t0=time.time()
    try:
        with contextlib.redirect_stdout(buf), warnings.catch_warnings():
            warnings.simplefilter("ignore")
            exec(compile(s,"mut","exec"),ns); ok=ns[fn]()
    except Exception as e:
        ok=f"EXC {type(e).__name__}: {e}"
    print(f"--- {label}  [{time.time()-t0:.0f}s]\n    -> returns {ok}\n    "+buf.getvalue().strip().replace("\n","\n    "),flush=True)
