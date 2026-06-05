#!/usr/bin/env python3
"""
1D Hou-Luo model (arXiv:2308.01528, arXiv:2604.01868):

    omega_t + u omega_x = theta_x   (+ nu omega_xx)
    theta_t + u theta_x = 0         (+ nu theta_xx)
    u_x = H[omega],   u(0) = 0   (u zero-mean on the periodic box)

Physical-space pseudospectral solver on a periodic box centered at the focus
x=0, domain [-pi, pi).  The model blows up in finite time at the origin.

Self-similar structure (smooth branch, 2308.01528):
    omega = (T-t)^{c_omega} Omega(xi),  theta = (T-t)^{c_theta} Theta(xi),
    xi = x/(T-t)^{c_l},   c_omega = -1,  c_theta = c_l - 2,  c_l in (2, 4.53).
    => max|omega| ~ (T-t)^{-1}.
Singular branch (2604.01868, Thm 2.3), triggered by the derivative-degeneracy
condition  omega_x(0)=0 and theta_xx(0)=0:  c_l -> 2 (endpoint), the rescaled
profile -> 1{X>1}/sqrt(X-1), far-field 1/sqrt(X).

Exponent extraction from a physical DNS (no a-priori exponents):
    c_omega : 1/max|omega| is linear in t with root at T   (c_omega=-1 test).
    c_l     : L(t) = max|omega| / max|omega_x| ~ (T-t)^{c_l}.
"""

import numpy as np
from numpy.fft import fft, ifft, fftfreq


class HouLuo:
    def __init__(self, N=4096, nu=0.0, L=2 * np.pi, dealias=True):
        self.N = int(N)
        self.nu = float(nu)
        self.L = float(L)
        self.dx = self.L / self.N
        # centered box [-L/2, L/2)
        self.x = -self.L / 2 + np.arange(self.N) * self.dx

        self.k = fftfreq(self.N, d=self.L / (2 * np.pi * self.N))  # integers for L=2pi
        self.ik = 1j * self.k
        self.k2 = self.k ** 2
        self.hilb = -1j * np.sign(self.k)
        self.hilb[0] = 0.0
        self.u_mult = np.zeros(self.N, dtype=complex)
        nz = self.k != 0
        self.u_mult[nz] = self.hilb[nz] / (1j * self.k[nz])  # u_hat = H(w)/(ik)

        if dealias:
            self.mask = (np.abs(self.k) <= self.N // 3).astype(float)
        else:
            self.mask = np.ones(self.N)

        self.omega = np.zeros(self.N)
        self.theta = np.zeros(self.N)

    def set_ic(self, omega0, theta0):
        self.omega = np.asarray(omega0, dtype=float).copy()
        self.theta = np.asarray(theta0, dtype=float).copy()
        return self

    def _dealias(self, f):
        return np.real(ifft(self.mask * fft(f)))

    def velocity(self, omega=None):
        if omega is None:
            omega = self.omega
        wh = fft(omega)
        u = np.real(ifft(self.u_mult * wh))
        u_x = np.real(ifft(self.hilb * wh))
        return u, u_x

    def rhs(self, omega, theta):
        wh = fft(omega)
        th = fft(theta)
        u = np.real(ifft(self.u_mult * wh))
        omega_x = np.real(ifft(self.ik * wh))
        theta_x = np.real(ifft(self.ik * th))

        r_omega = self._dealias(-u * omega_x) + theta_x
        r_theta = self._dealias(-u * theta_x)
        if self.nu > 0.0:
            r_omega = r_omega + self.nu * np.real(ifft(-self.k2 * wh))
            r_theta = r_theta + self.nu * np.real(ifft(-self.k2 * th))
        return r_omega, r_theta

    def step_rk4(self, dt):
        w, t = self.omega, self.theta
        k1w, k1t = self.rhs(w, t)
        k2w, k2t = self.rhs(w + 0.5 * dt * k1w, t + 0.5 * dt * k1t)
        k3w, k3t = self.rhs(w + 0.5 * dt * k2w, t + 0.5 * dt * k2t)
        k4w, k4t = self.rhs(w + dt * k3w, t + dt * k3t)
        self.omega = w + (dt / 6.0) * (k1w + 2 * k2w + 2 * k3w + k4w)
        self.theta = t + (dt / 6.0) * (k1t + 2 * k2t + 2 * k3t + k4t)

    def auto_dt(self, cfl=0.1, dt_max=1e-2):
        u, u_x = self.velocity()
        umax = np.max(np.abs(u)) + 1e-30
        dt = cfl * self.dx / umax            # advective CFL
        dt = min(dt, cfl / (np.max(np.abs(u_x)) + 1e-30))  # stretching rate
        if self.nu > 0.0:
            dt = min(dt, 0.4 * self.dx ** 2 / (2 * self.nu))
        return min(dt, dt_max)

    def tail_fraction(self):
        """Resolution alarm = energy in the UPPER RETAINED band |k| in [N/6, N/3].
        (Modes |k|>N/3 are zeroed by the 2/3 dealias mask, so they are useless as
        an alarm -- they read ~0 by construction.  When the structure sharpens
        toward grid scale, the retained band near the cutoff fills first.)"""
        p = np.abs(fft(self.omega)) ** 2 + np.abs(fft(self.theta)) ** 2
        band = (np.abs(self.k) >= self.N // 6) & (np.abs(self.k) <= self.N // 3)
        return float(np.sum(p[band]) / (np.sum(p) + 1e-30))

    def length_scale(self):
        """L = max|omega| / max|omega_x| ~ (T-t)^{c_l}."""
        wh = fft(self.omega)
        omega_x = np.real(ifft(self.ik * wh))
        return float(np.max(np.abs(self.omega)) / (np.max(np.abs(omega_x)) + 1e-30))

    def diagnostics(self):
        wmax = float(np.max(np.abs(self.omega)))
        loc = int(np.argmax(np.abs(self.omega)))
        return {
            "omega_max": wmax,
            "omega_argmax_x": float(self.x[loc]),
            "theta_max": float(np.max(np.abs(self.theta))),
            "length_scale": self.length_scale(),
            "tail_fraction": self.tail_fraction(),
        }


def run_physical(solver, T_max, cfl=0.1, dt_max=1e-2, sample_every=25,
                 blowup_thresh=1e7, tail_thresh=0.05, verbose=False):
    """Integrate until T_max / blowup / under-resolution.  Returns history."""
    t = 0.0
    step = 0
    hist = []
    status = "reached_T"
    while t < T_max:
        d = solver.diagnostics()
        if d["omega_max"] > blowup_thresh:
            status = "blowup"
            d["time"] = t
            hist.append(d)
            break
        if d["omega_max"] > 5.0 and d["tail_fraction"] > tail_thresh:
            status = "under_resolved"
            d["time"] = t
            hist.append(d)
            break
        if step % sample_every == 0:
            d["time"] = t
            hist.append(d)
            if verbose and step % (sample_every * 20) == 0:
                print(f"    t={t:.5f}  |w|_max={d['omega_max']:.4e}  "
                      f"L={d['length_scale']:.4e}  tail={d['tail_fraction']:.1e}")
        dt = solver.auto_dt(cfl=cfl, dt_max=dt_max)
        dt = min(dt, T_max - t)
        solver.step_rk4(dt)
        t += dt
        step += 1
    return {"status": status, "t_final": t, "n_steps": step, "history": hist}


def fit_self_similar(history, tail_frac=0.4):
    """Estimate (T, c_omega, c_l) from the tail of a blowing-up run.

    c_omega from 1/max|omega| ~ (1/A)(T-t): linear fit of 1/M vs t, root = T,
        and a log-log check that M ~ (T-t)^{c_omega}.
    c_l from log L vs log(T-t), L = length_scale.
    """
    t = np.array([h["time"] for h in history])
    M = np.array([h["omega_max"] for h in history])
    Lsc = np.array([h["length_scale"] for h in history])
    n = len(t)
    i0 = int(n * (1 - tail_frac))
    tt, MM, LL = t[i0:], M[i0:], Lsc[i0:]

    # T from linear fit of 1/M vs t (valid iff c_omega = -1).
    inv = 1.0 / MM
    A = np.polyfit(tt, inv, 1)         # inv ~ a*t + b ; root T = -b/a
    T_lin = -A[1] / A[0]

    # log-log slopes against (T - t) using T_lin.
    tau = T_lin - tt
    good = tau > 0
    c_omega = np.polyfit(np.log(tau[good]), np.log(MM[good]), 1)[0]
    c_l = np.polyfit(np.log(tau[good]), np.log(LL[good]), 1)[0]
    return {"T": float(T_lin), "c_omega": float(c_omega), "c_l": float(c_l),
            "n_fit": int(np.sum(good))}
