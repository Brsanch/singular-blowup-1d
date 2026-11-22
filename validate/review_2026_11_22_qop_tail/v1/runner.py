"""argv: script tail(0/1/default). Runs repo script with build_blocks_z's edge_tail default forced (0/1) or untouched."""
import sys, os, runpy
repo = "/Volumes/4TB SD/ClaudeCode/singular-blowup-1d"
sys.path.insert(0, os.path.join(repo, "src")); sys.path.insert(0, os.path.join(repo, "validate"))
import stability_edge_z as m
mode = sys.argv[2]
if mode in ("0", "1"):
    orig = m.build_blocks_z
    def patched(*a, **k):
        k.setdefault("edge_tail", mode == "1"); return orig(*a, **k)
    m.build_blocks_z = patched
script = os.path.join(repo, sys.argv[1]); sys.argv = [script]; os.chdir("/private/tmp/qop_review/cwd")
import matplotlib; matplotlib.use("Agg")
runpy.run_path(script, run_name="__main__")
