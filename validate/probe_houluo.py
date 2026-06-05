#!/usr/bin/env python3
"""
Exploratory probe: does the Hou-Luo solver blow up, and what is c_omega?

Streams a flushed trajectory so growth/decay is visible live.  Stops at modest
amplification (enough to fit exponents) or on decay.  Parity: omega odd,
theta even about the focus x=0; box [-pi, pi).
"""
import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from houluo import HouLuo, fit_self_similar  # noqa: E402


def make_ic(name, x, amp):
    if name == "sin_cos":
        return amp * np.sin(x), amp * np.cos(x)
    if name == "sin_negcos":
        return amp * np.sin(x), -amp * np.cos(x)
    if name == "sin_1mcos":
        return amp * np.sin(x), amp * (1.0 - np.cos(x))
    raise ValueError(name)


def run_stream(s, amp_stop=2.0e3, max_steps=200000, cfl=0.08, dt_max=2e-2,
               print_decades=True):
    """Manual loop with flushed streaming. Stops at amp_stop or decay."""
    t = 0.0
    hist = []
    w0max = float(np.max(np.abs(s.omega)))
    next_print = w0max * 2.0
    decay_count = 0
    status = "max_steps"
    for step in range(max_steps):
        d = s.diagnostics()
        d["time"] = t
        wmax = d["omega_max"]
        if step % 25 == 0:
            hist.append(d)
        if wmax > amp_stop:
            status = "amp_stop"
            hist.append(d)
            break
        if wmax > 5.0 and d["tail_fraction"] > 0.05:
            status = "under_resolved"
            hist.append(d)
            break
        if print_decades and wmax >= next_print:
            print(f"      t={t:.5f}  |w|max={wmax:.4e}  L={d['length_scale']:.3e}  "
                  f"x*={d['omega_argmax_x']:+.4f}  tail={d['tail_fraction']:.1e}",
                  flush=True)
            next_print *= 2.0
        # decay detection
        if wmax < 0.4 * w0max:
            decay_count += 1
            if decay_count > 2000:
                status = "decayed"
                break
        else:
            decay_count = 0
        dt = s.auto_dt(cfl=cfl, dt_max=dt_max)
        s.step_rk4(dt)
        t += dt
    return {"status": status, "t_final": t, "history": hist, "w0max": w0max,
            "wmax_final": float(np.max(np.abs(s.omega)))}


def main():
    N = 2048
    cases = [("sin_cos", 1.0), ("sin_negcos", 1.0), ("sin_1mcos", 1.0)]
    print("=" * 72, flush=True)
    print("  HOU-LUO PHYSICAL-DNS PROBE  (smooth/non-degenerate ICs, N=%d)" % N,
          flush=True)
    print("=" * 72, flush=True)
    for name, amp in cases:
        s = HouLuo(N=N, nu=0.0)
        w0, t0 = make_ic(name, s.x, amp)
        s.set_ic(w0, t0)
        print(f"\n  IC={name:10s} amp={amp}  |w|max0={s.diagnostics()['omega_max']:.3f}",
              flush=True)
        res = run_stream(s)
        print(f"    -> status={res['status']}  t_final={res['t_final']:.4f}  "
              f"|w|max {res['w0max']:.3f} -> {res['wmax_final']:.4e}", flush=True)
        h = res["history"]
        if res["status"] in ("amp_stop", "under_resolved") and len(h) > 12:
            fit = fit_self_similar(h, tail_frac=0.4)
            print(f"    FIT  T={fit['T']:.5f}  c_omega={fit['c_omega']:+.4f} "
                  f"(expect -1)  c_l={fit['c_l']:+.4f} (expect (2,4.53))  "
                  f"n_fit={fit['n_fit']}", flush=True)
    print("\n" + "=" * 72, flush=True)


if __name__ == "__main__":
    main()
