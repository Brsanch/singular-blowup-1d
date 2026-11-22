import sys, importlib.util, numpy as np
REPO="/Volumes/4TB SD/ClaudeCode/singular-blowup-1d"
def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path); m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m
NEW=load(REPO+"/src/stability_edge_z.py","sez_new")
AD4=load("/private/tmp/hdiag_review/old_ad4/stability_edge_z.py","sez_ad4")
PRE=load("/private/tmp/hdiag_review/old_pre/stability_edge_z.py","sez_pre")
