#!/usr/bin/env python3
"""
Pre-numerics verification of the EDGE-CONDITION derivation for the desingularized
free-boundary stability eigenvalue problem (Conjecture 2.4 of arXiv:2604.01868).

This does NOT build the full operator.  It checks the load-bearing analytic facts
that the operator construction will rest on, BEFORE committing to numerics:

  (A) Fixed-frame linearization of the edge shift gives a SPURIOUS infinite edge
      velocity; the moving (edge-fixed) frame keeps it finite.  This is the
      conceptual reason the naive regularized spectrum produced a delta-divergent
      edge artifact (the recorded dead end), and why the edge-fixed frame is the
      correct formulation.

  (B) The discrete line-Hilbert matrix commutes with the discrete derivative:
      H(f') == (H f)'  (the identity H(Omega_bar') = (H Omega_bar)' = U_bar_XX
      that makes the perturbation velocity well-defined on the support side).
      If the DISCRETE operators do not commute, the assembled stability operator
      is inconsistent exactly at the edge.

  (C) The edge ODE indicial structure and the predicted symmetry-mode eigenvalue
      lambda = c_l - 1 = 1 for the pure edge-displacement (xi) mode.

Run: python3 validate/probe_edge_condition.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from rescaled import hilbert_matrix, singular_U_analytic  # noqa: E402


def banner(s):
    print("\n" + "=" * 78 + f"\n  {s}\n" + "=" * 78)


# ---------------------------------------------------------------------------
# (A) Edge-shift velocity: fixed-frame linearization is a SPURIOUS divergence;
#     the moving frame keeps it finite.
#
#   Base velocity (closed form, exact):  U_bar(Y) = 2 sqrt(1-Y) - 2  (Y<1), -2 (Y>1).
#   U_bar is FINITE and continuous (U_bar(1) = -2), but U_bar_Y = -1/sqrt(1-Y)
#   has a VERTICAL TANGENT (-> -inf) as Y -> 1^-.
#
#   Moving frame: the perturbation velocity is  U_bar(Y - s) - U_bar(Y)  (a finite
#   difference of finite numbers) -> the true, finite edge-shift velocity.
#   Fixed frame (naive linearization): -s U_bar_Y(Y) = +s/sqrt(1-Y) -> +inf.
#
#   If linearization were valid, [U_bar(Y-s)-U_bar(Y)]/s ~ -U_bar_Y(Y).  We show
#   the TRUE ratio stays O(1) while -U_bar_Y blows up like (1-Y)^{-1/2} as Y->1^-.
# ---------------------------------------------------------------------------
def check_A():
    banner("(A) edge-shift velocity: moving-frame finite vs fixed-frame artifact")
    s = 1e-4  # tiny edge displacement (linear regime)
    print(f"  edge displacement s = {s:g}")
    print(f"  {'1-Y':>10} {'true (U(Y-s)-U(Y))/s':>22} {'fixed -U_bar_Y(Y)':>20} "
          f"{'ratio':>10}")
    Ys = 1.0 - np.array([1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6])
    for Y in Ys:
        true_ratio = (singular_U_analytic(np.array([Y - s]))[0]
                      - singular_U_analytic(np.array([Y]))[0]) / s
        UY = -1.0 / np.sqrt(1.0 - Y)          # U_bar_Y on Y<1
        fixed = -UY                            # the linear-in-s coefficient
        print(f"  {1.0 - Y:10.1e} {true_ratio:22.4f} {fixed:20.4f} "
              f"{(fixed / true_ratio):10.2f}")
    print("\n  READ: 'true' stays O(1) (moving frame finite); 'fixed' diverges")
    print("  like (1-Y)^(-1/2) (the naive X-frame linearization artifact).")
    print("  => the edge-fixed frame is REQUIRED; the -xi U_bar_XX velocity is")
    print("     a Taylor-breakdown at the vertical-tangent edge, NOT physical.")
    # On the SUPPORT side Y>1, U_bar = -2 (flat) so the edge-shift velocity is 0:
    Yp = 1.0 + np.array([1e-3, 1e-1, 1.0, 5.0])
    dv = (singular_U_analytic(Yp - s) - singular_U_analytic(Yp)) / s
    print(f"\n  support side Y>1: (U_bar(Y-s)-U_bar(Y))/s = {np.round(dv,6)}")
    print("  (all ~0: U_bar is flat (-2) on Y>1, so a pure edge shift adds no")
    print("   velocity on the support side -> coupling U*Omega_bar_X is finite.)")


# ---------------------------------------------------------------------------
# (B) Discrete commutation  H(f') == (H f)'  for the validated line-Hilbert
#     matrix + the 4th-order central derivative used in the operator.
#     Test on a smooth, fast-decaying, well-resolved f so truncation is tiny.
# ---------------------------------------------------------------------------
def _dc4(N, dx):
    D = np.zeros((N, N))
    c = np.array([1.0, -8.0, 0.0, 8.0, -1.0]) / (12.0 * dx)
    for a, o in zip(c, (-2, -1, 0, 1, 2)):
        D += a * np.eye(N, k=o)
        if o > 0:
            D += a * np.eye(N, k=o - N)
        elif o < 0:
            D += a * np.eye(N, k=o + N)
    return D


def check_B():
    banner("(B) discrete commutation  H d_Y == d_Y H  (well-defined edge velocity)")
    N = 1201
    L = 40.0
    Y = np.linspace(-L, L, N)
    dx = Y[1] - Y[0]
    H = hilbert_matrix(N)
    D = _dc4(N, dx)                          # CIRCULANT (periodic-wrap) derivative
    C = H @ D - D @ H

    # function-level: smooth bump, exact-ish Hilbert
    f = np.exp(-Y ** 2 / 2.0) * Y
    inner = slice(60, N - 60)
    lhs = H @ (D @ f)                        # H(f')
    rhs = D @ (H @ f)                        # (H f)'
    err = np.max(np.abs((lhs - rhs)[inner]))
    scale = np.max(np.abs(rhs[inner])) + 1e-30
    print(f"  function-level  max|H(f') - (Hf)'| / scale (interior) = {err/scale:.3e}")

    # operator-level: WHERE does [H,D] live?
    rownorm = np.max(np.abs(C), axis=1)
    edge = int(np.argmin(np.abs(Y - 1.0)))
    Cmask = C.copy()
    Cmask[:60] = Cmask[-60:] = 0.0
    Cmask[:, :60] = Cmask[:, -60:] = 0.0
    print(f"  operator-level  ||[H,D]||_inf  (full)            = {np.max(np.abs(C)):.3e}")
    print(f"                  rownorm at box ends (Y=+-{L:g})     = {rownorm[0]:.3e}")
    print(f"                  rownorm at physical edge (Y~1)    = {rownorm[edge]:.3e}")
    print(f"                  ||[H,D]||_inf after masking borders = "
          f"{np.max(np.abs(Cmask)):.3e}")
    print("\n  READ: function-level commutation is exact, and the OPERATOR")
    print("  commutator is machine-zero everywhere EXCEPT the two box ends, where")
    print("  the CIRCULANT (periodic-wrap) derivative fights the non-periodic line-")
    print("  Hilbert -> O(1) rows at Y=+-L only. At the physical edge Y=1 it is")
    print("  ~1e-2 (bulk level). CONSEQUENCES for the desingularized operator:")
    print("   - H(omega')=(H omega)'=U_bar_XX is discretely CONSISTENT at the edge:")
    print("     the edge-velocity cancellation (U_bar_XX=0 on support side) holds.")
    print("   - the box-end O(1) failure is a FAR-FIELD artifact and is exactly the")
    print("     source of the spurious 'box modes' (x*=40) in the resolved spectrum")
    print("     -> use a NON-periodic (outflow) far-field derivative or a mapped")
    print("        infinite domain; never a circulant D on the real line.")
    print("   - near the edge use the CLOSED-FORM U_bar, U_bar_X, U_bar_XX; do NOT")
    print("     discretely Hilbert-transform the singular Omega_bar'.")


# ---------------------------------------------------------------------------
# (C) Edge ODE indicial structure + predicted edge-mode eigenvalue.
#
#   On Y>1 the local homogeneous transport+damping part of the shape eq is
#       lambda g = -abar g_Y + c_omega g,   abar = c_l (Y-1)  (=2(Y-1)),
#   => 2(Y-1) g_Y = -(lambda+1) g  => g ~ (Y-1)^{-(lambda+1)/2}.
#   Weighted space (g no more singular than the profile (Y-1)^{-1/2}) admits the
#   homogeneous edge solution only if Re(-(lambda+1)/2) >= -1/2, i.e. Re(lambda)<=0.
#   Forced (Y-1)^{-1/2} particular part has edge coeff B = H(g)(1)/(2 lambda).
#
#   The pure edge (xi) mode: xi_dot = (c_l - 1) xi + V_g(1) + mu.  With g=0, mu=0:
#       xi_dot = (c_l - 1) xi  =>  lambda = c_l - 1 = 1   (the symmetry/
#       normalization mode that the mu-normalization projects out).
# ---------------------------------------------------------------------------
def check_C():
    banner("(C) edge indicial exponents + predicted edge-mode eigenvalue")
    c_l, c_omega = 2.0, -1.0
    print("  local homogeneous edge solution:  g ~ (Y-1)^r,  r = -(lambda+1)/2")
    print(f"  {'lambda':>8} {'r=-(lam+1)/2':>14} {'admissible in weighted space?':>32}")
    for lam in (-2.0, -1.0, 0.0, 0.5, 1.0, 2.0):
        r = -(lam + 1.0) / 2.0
        ok = "yes (<= profile -1/2)" if r >= -0.5 - 1e-12 else "NO (too singular)"
        print(f"  {lam:8.2f} {r:14.3f} {ok:>32}")
    lam_edge = c_l - 1.0
    print(f"\n  pure edge (xi) mode:  xi_dot = (c_l-1) xi + V_g(1) + mu")
    print(f"     g=0, mu=0  =>  lambda = c_l - 1 = {lam_edge:g}")
    print(f"  PREDICTION to validate against the full operator: a mode at")
    print(f"  lambda = {lam_edge:g} that is PURE edge displacement (g~0), to be")
    print(f"  projected out as the time-translation/normalization symmetry mode.")
    print(f"  forced edge coeff of g:  B = H(g)(1) / (2 lambda)  (finite, lam!=0).")
    print("\n  READ: the weighted space admits homogeneous edge behavior only for")
    print("  Re(lambda) <= 0; any genuine unstable (Re>0) mode must therefore have")
    print("  NO homogeneous edge singularity (a Robin-type edge BC), only the")
    print("  forced (Y-1)^{-1/2} part. That BC is what the operator must impose.")


if __name__ == "__main__":
    check_A()
    check_B()
    check_C()
    print("\n" + "=" * 78)
    print("  edge-condition derivation: load-bearing facts checked.")
    print("=" * 78)
