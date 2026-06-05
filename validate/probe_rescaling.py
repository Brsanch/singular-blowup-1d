#!/usr/bin/env python3
"""
Probe: does the dynamic-rescaling solver converge to a smooth self-similar
profile with c_omega -> -1 for a prescribed c_l?

Validation target (arXiv:2308.01528): smooth profiles exist for c_l in (2,4.53),
all with c_omega = -1.  So c_omega should converge to -1 for c_l=3 (inside) and
should fail / not settle at -1 for c_l outside the interval.
"""
import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from dynamic_rescaling import DynamicRescaling, run  # noqa: E402


def smooth_ic(X):
    Omega0 = X * np.exp(-X ** 2 / 2.0)     # odd, Omega_X(0)=1
    Theta0 = np.exp(-X ** 2 / 2.0)         # even, Theta(0)=1
    return Omega0, Theta0


def main():
    N, P = 1025, 60.0
    print("=" * 74)
    print("  DYNAMIC-RESCALING PROBE  (smooth branch; target c_omega->-1)")
    print(f"  N={N}, P={P}")
    print("=" * 74)
    for c_l in (3.0,):
        print(f"\n  --- prescribed c_l = {c_l} ---", flush=True)
        s = DynamicRescaling(N=N, P=P, c_l=c_l)
        s.set_ic(*smooth_ic(s.X))
        hist = run(s, tau_max=12.0, cfl=0.3, dt_max=5e-3, log_every=1500)
        if hist:
            last = hist[-1]
            print(f"\n  RESULT c_l={c_l}: c_omega={last['c_omega']:+.5f} "
                  f"(target -1), residual={last['residual']:.3e}", flush=True)


if __name__ == "__main__":
    main()
