"""Review script (2026-11-22, post-f5106e9 review). Generic mutation runner: python3 mut.py <gate.py> <check_fn> <label> <old>|||<new> [<old>|||<new> ...]"""
import sys, io, contextlib, traceback
path, fn, label = sys.argv[1], sys.argv[2], sys.argv[3]
src = open(path).read()
for spec in sys.argv[4:]:
    old, new = spec.split("|||")
    n = src.count(old)
    if n == 0:
        print(f"[{label}] MUTATION NOT APPLIED: {old!r}"); sys.exit(2)
    src = src.replace(old, new)
ns = {"__name__": "mutant", "__file__": path}
exec(compile(src, path, "exec"), ns)
buf = io.StringIO()
try:
    with contextlib.redirect_stdout(buf):
        r = ns[fn]()
    out = buf.getvalue().strip().splitlines()
    print(f"[{label}] returned PASS={r} :: {out[-1][:260] if out else ''}")
except Exception as e:
    print(f"[{label}] EXCEPTION {type(e).__name__}: {str(e)[:200]}")
