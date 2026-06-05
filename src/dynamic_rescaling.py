#!/usr/bin/env python3
"""
Dynamic-rescaling (self-similar coordinate) solver for the 1D Hou-Luo model,
eq. (2.4) of arXiv:2604.01868:

    Omega_tau + (U + c_l X) Omega_X = c_omega Omega + Theta_X
    Theta_tau + (U + c_l X) Theta_X = (c_l + 2 c_omega) Theta
    U_X = H(Omega),  U(0) = 0      (Hilbert transform on R)

Steady states are the self-similar profiles.  Smooth (non-degenerate) family
exists for c_l in (2, 4.53), all with c_omega = -1 (arXiv:2308.01528); the
endpoint c_l=2 is the singular profile Omega_bar = 1{X>1}/sqrt(X-1).

NUMERICS (chosen to respect the structure):
  * Domain R approximated by a large box [-P, P], equispaced (for the line
    Hilbert transform).  Profiles decay only algebraically (~X^{-1/c_l}), so the
    box must be large; convergence is checked by varying P.
  * The transport term (U + c_l X) d_X is advected with a 3rd-order UPWIND
    scheme.  The characteristic speed a = U + c_l X ~ c_l X is OUTGOING at both
    ends (a>0 for X>>0, a<0 for X<<0), so upwind points inward -> stable
    outflow, no inflow BC needed.  (A periodic/Fourier method is INVALID here:
    c_l X is not periodic.)
  * Source terms (+Theta_X in the Omega-eq) use 4th-order CENTRAL differences.
  * U_X = H(Omega) via the validated line Hilbert (src.rescaled.hilbert_matrix);
    U by cumulative integration anchored at U(0)=0.
  * Normalization: prescribe c_l; hold Omega_X(0)=p fixed.  d_X of the Omega-eq
    at X=0 gives c_omega = h0 + c_l - Theta_XX(0)/p,  h0 = H(Omega)(0).
    At a converged smooth profile c_omega -> -1 (the validation target).
"""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.dirname(__file__))
from rescaled import hilbert_matrix  # noqa: E402


def _roll(a, s):
    return np.roll(a, s)


def d_central4(f, dx):
    return (-_roll(f, -2) + 8 * _roll(f, -1) - 8 * _roll(f, 1) + _roll(f, 2)) / (12 * dx)


def d2_central(f, dx):
    return (_roll(f, -1) - 2 * f + _roll(f, 1)) / dx ** 2


def d_upwind3(f, a, dx):
    """3rd-order upwind-biased d_X f, branch by sign of speed a."""
    # a>0: bias to the left (points i-2,i-1,i,i+1)
    dpos = (2 * _roll(f, -1) + 3 * f - 6 * _roll(f, 1) + _roll(f, 2)) / (6 * dx)
    # a<0: bias to the right (points i-1,i,i+1,i+2)
    dneg = (-_roll(f, -2) + 6 * _roll(f, -1) - 3 * f - 2 * _roll(f, 1)) / (6 * dx)
    return np.where(a >= 0.0, dpos, dneg)


class DynamicRescaling:
    def __init__(self, N=1025, P=60.0, c_l=3.0, c_omega_fixed=None):
        assert N % 2 == 1, "use odd N so X=0 is a grid node"
        self.N = N
        self.P = P
        self.c_l = float(c_l)
        # If set, use this fixed c_omega instead of the X=0 constraint (the smooth
        # family has c_omega=-1 for all c_l in (2,4.53)).
        self.c_omega_fixed = c_omega_fixed
        self.X = np.linspace(-P, P, N)
        self.dx = self.X[1] - self.X[0]
        self.c = N // 2                      # index of X=0
        assert abs(self.X[self.c]) < 1e-12
        self.M = hilbert_matrix(N)           # line Hilbert (dense)
        self.Omega = np.zeros(N)
        self.Theta = np.zeros(N)
        self.p0 = None                       # fixed Omega_X(0)
        self.c_omega = -1.0

    def set_ic(self, Omega0, Theta0):
        self.Omega = np.asarray(Omega0, float).copy()
        self.Theta = np.asarray(Theta0, float).copy()
        self.p0 = float(d_central4(self.Omega, self.dx)[self.c])
        return self

    def velocity(self, Omega):
        Ux = self.M @ Omega                  # H(Omega)
        # U(X) = int_0^X Ux dX', U(0)=0  (cumulative trapezoid, anchored at c)
        U = np.concatenate(([0.0], np.cumsum(0.5 * (Ux[1:] + Ux[:-1]) * self.dx)))
        U = U - U[self.c]
        return U, Ux

    def c_omega_from_constraint(self, Omega, Theta, Ux):
        h0 = Ux[self.c]                       # H(Omega)(0)
        q = d2_central(Theta, self.dx)[self.c]   # Theta_XX(0)
        return h0 + self.c_l - q / self.p0

    def rhs(self, Omega, Theta):
        U, Ux = self.velocity(Omega)
        a = U + self.c_l * self.X
        if self.c_omega_fixed is not None:
            c_omega = self.c_omega_fixed
        else:
            c_omega = self.c_omega_from_constraint(Omega, Theta, Ux)
        Om_X = d_upwind3(Omega, a, self.dx)
        Th_X_adv = d_upwind3(Theta, a, self.dx)
        Th_X_src = d_central4(Theta, self.dx)
        rO = -a * Om_X + c_omega * Omega + Th_X_src
        rT = -a * Th_X_adv + (self.c_l + 2 * c_omega) * Theta
        return rO, rT, c_omega

    def step_rk4(self, dt):
        w, t = self.Omega, self.Theta
        k1w, k1t, c1 = self.rhs(w, t)
        k2w, k2t, _ = self.rhs(w + 0.5 * dt * k1w, t + 0.5 * dt * k1t)
        k3w, k3t, _ = self.rhs(w + 0.5 * dt * k2w, t + 0.5 * dt * k2t)
        k4w, k4t, _ = self.rhs(w + dt * k3w, t + dt * k3t)
        self.Omega = w + dt / 6 * (k1w + 2 * k2w + 2 * k3w + k4w)
        self.Theta = t + dt / 6 * (k1t + 2 * k2t + 2 * k3t + k4t)
        self.c_omega = c1
        return c1

    def auto_dt(self, cfl=0.3, dt_max=5e-3):
        U, _ = self.velocity(self.Omega)
        amax = np.max(np.abs(U + self.c_l * self.X)) + 1e-30
        return min(cfl * self.dx / amax, dt_max)

    def residual_norm(self):
        rO, rT, _ = self.rhs(self.Omega, self.Theta)
        # focus region only (the far field is slow to settle / box-limited)
        core = np.abs(self.X) < 8.0
        return float(np.max(np.abs(rO[core])) + np.max(np.abs(rT[core])))

    def omega_x0(self):
        return float(d_central4(self.Omega, self.dx)[self.c])


def run(solver, tau_max, cfl=0.3, dt_max=5e-3, log_every=2000, verbose=True):
    tau = 0.0
    step = 0
    hist = []
    while tau < tau_max:
        dt = solver.auto_dt(cfl=cfl, dt_max=dt_max)
        c_om = solver.step_rk4(dt)
        tau += dt
        step += 1
        if step % log_every == 0:
            res = solver.residual_norm()
            rec = {"tau": tau, "c_omega": c_om, "residual": res,
                   "omega_x0": solver.omega_x0(),
                   "max_omega": float(np.max(np.abs(solver.Omega)))}
            hist.append(rec)
            if verbose:
                print(f"    tau={tau:7.3f} c_omega={c_om:+.5f} resid={res:.3e} "
                      f"Om_X(0)={rec['omega_x0']:+.4f} max|Om|={rec['max_omega']:.3f}",
                      flush=True)
        if not np.all(np.isfinite(solver.Omega)):
            print("    BLEW UP (non-finite)", flush=True)
            break
    return hist
