import re,sys,collections
d=collections.OrderedDict()
for l in open(sys.argv[1]):
    m=re.match(r"N=\s*(\d+) zmin=\s*(-?\d+) zmax=\s*(\d+) bc=(\d) g=([\d.]+)\s+repo_tail\(Q_op/2\)=(\S+) repo_tail\(forcing\)=(\S+) EXACTKERNEL=(\S+) truncated=\s*(\S+) \| argmax\|v\| z: repo (\S+) exact (\S+)",l)
    if m:
        k=m.group(1,2,3,4); d.setdefault(k,[]).append(m.groups()[4:])
for k,v in d.items():
    print("N=%s zmin=%s zmax=%s bc=%s"%k)
    print("   g:    "+" ".join(f"{x[0]:>9s}" for x in v))
    print("   repo: "+" ".join(f"{float(x[1]):+9.2e}" for x in v))
    print("   exact:"+" ".join(f"{float(x[3]):+9.2e}" for x in v))
    print("   trunc:"+" ".join(f"{float(x[4]):+9.2f}" for x in v))
    print("   zrepo:"+" ".join(f"{x[5]:>9s}" for x in v))
    print("   zexct:"+" ".join(f"{x[6]:>9s}" for x in v))
    print("   max|Q/2 - forcing| FoV diff: %.1e"%max(abs(float(x[1])-float(x[2])) for x in v))
