import sys, runpy, importlib.util
spec=importlib.util.spec_from_file_location("stability_edge_z","/private/tmp/hdiag_review/old_ad4/stability_edge_z.py")
m=importlib.util.module_from_spec(spec); sys.modules["stability_edge_z"]=m; spec.loader.exec_module(m)
script=sys.argv[1]; sys.argv=[script]
import os; sys.path.insert(0,os.path.dirname(script))
runpy.run_path(script,run_name="__main__")
