#!/usr/bin/env python3
"""
Generalized Constantin-Lax-Majda (gCLM) / De Gregorio pseudospectral solver.

    omega_t + a * u * omega_x = u_x * omega + nu * omega_xx
    u_x = H[omega],     u zero-mean (periodic)

on the circle x in [0, 2*pi).  H is the Hilbert transform with the standard
convention  H[e^{i k x}] = -i * sign(k) * e^{i k x}  (so H[sin] = -cos,
H[cos] = sin).

Special cases of the advection parameter a:
    a = 0   ->  Constantin-Lax-Majda (CLM, 1985):  omega_t = omega * H[omega]
               (no advection; exactly solvable -- see clm_exact.py)
    a = 1   ->  De Gregorio (1990):  omega_t + u omega_x = u_x omega
    general a -> Okamoto-Sakajo-Wunsch (2008) gCLM family.

Sign convention for the stretching term is  +u_x * omega  (vortex stretching),
matching the OSW/De Gregorio literature and the CLM exact solution.

This is the shared engine for the singular-self-similar-blowup program:
CLM/De Gregorio/gCLM all run through this class; the Hou-Luo coupled model
lives in houluo.py and reuses the same spectral primitives.

Design notes
------------
- All spectral.  Quadratic nonlinearities (u*omega_x, u_x*omega) are formed in
  physical space then 2/3-dealiased, which is exact for quadratic products.
- RK4 in time with adaptive dt driven by the stretching rate max|u_x| and the
  advective CFL.  Viscosity (nu>0) handled inside the RHS (explicit); for the
  inviscid blowup studies nu=0 and the limiter is purely the stretching rate.
- Everything is float64.  N is a power of two for FFT speed; 1D so N up to ~2^22
  is cheap on this machine.
"""

import numpy as np
from numpy.fft import rfft, irfft, fft, ifft, fftfreq


class GCLM:
    def __init__(self, N=2048, a=0.0, nu=0.0, L=2 * np.pi, dealias=True):
        self.N = int(N)
        self.a = float(a)
        self.nu = float(nu)
        self.L = float(L)
        self.dx = self.L / self.N
        self.x = np.arange(self.N) * self.dx  # [0, L)

        # Integer wavenumbers for L = 2*pi.  fftfreq(N, d=L/(2 pi N)) = integers.
        self.k = fftfreq(self.N, d=self.L / (2 * np.pi * self.N))
        self.ik = 1j * self.k
        self.k2 = self.k ** 2

        # Hilbert multiplier  H_hat(k) = -i sign(k), with H(0)=0.
        self.hilb = -1j * np.sign(self.k)
        self.hilb[0] = 0.0

        # 1/(ik) for recovering u from u_x = H[omega]  (u zero-mean).
        self.inv_ik = np.zeros(self.N, dtype=complex)
        nz = self.k != 0
        self.inv_ik[nz] = 1.0 / (1j * self.k[nz])

        # u_hat = hilb * omega_hat / (i k)  ->  precompute the combined multiplier.
        self.u_mult = self.hilb * self.inv_ik  # = -i sign(k) / (i k) = -sign(k)/k

        # 2/3 dealias mask.
        if dealias:
            self.mask = np.ones(self.N)
            cutoff = self.N // 3
            self.mask[np.abs(self.k) > cutoff] = 0.0
        else:
            self.mask = np.ones(self.N)

        self.omega = np.zeros(self.N)

    # ---- field setters -------------------------------------------------
    def set_omega(self, omega0):
        self.omega = np.asarray(omega0, dtype=float).copy()
        return self

    # ---- spectral primitives ------------------------------------------
    def hilbert(self, f):
        return np.real(ifft(self.hilb * fft(f)))

    def deriv(self, f):
        return np.real(ifft(self.ik * fft(f)))

    def velocity(self, omega=None):
        """Return (u, u_x) with u_x = H[omega], u zero-mean."""
        if omega is None:
            omega = self.omega
        wh = fft(omega)
        u = np.real(ifft(self.u_mult * wh))
        u_x = np.real(ifft(self.hilb * wh))
        return u, u_x

    # ---- dynamics ------------------------------------------------------
    def rhs(self, omega):
        wh = fft(omega)
        u = np.real(ifft(self.u_mult * wh))
        u_x = np.real(ifft(self.hilb * wh))
        omega_x = np.real(ifft(self.ik * wh))

        r = -self.a * (u * omega_x) + (u_x * omega)
        # 2/3-dealias the quadratic products.
        r = np.real(ifft(self.mask * fft(r)))

        if self.nu > 0.0:
            r = r + self.nu * np.real(ifft(-self.k2 * wh))
        return r

    def step_rk4(self, dt):
        w = self.omega
        k1 = self.rhs(w)
        k2 = self.rhs(w + 0.5 * dt * k1)
        k3 = self.rhs(w + 0.5 * dt * k2)
        k4 = self.rhs(w + dt * k3)
        self.omega = w + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

    def auto_dt(self, cfl=0.1, dt_max=1e-2):
        """dt from stretching rate max|u_x| and advective CFL max|u|/dx."""
        u, u_x = self.velocity()
        rate = np.max(np.abs(u_x))
        dt = cfl / (rate + 1e-30)
        if self.a != 0.0:
            umax = np.max(np.abs(u))
            dt_adv = cfl * self.dx / (umax + 1e-30)
            dt = min(dt, dt_adv)
        if self.nu > 0.0:
            dt = min(dt, 0.4 * self.dx ** 2 / (2 * self.nu))
        return min(dt, dt_max)

    # ---- diagnostics ---------------------------------------------------
    def omega_max(self):
        return float(np.max(np.abs(self.omega)))

    def tail_fraction(self):
        """Fraction of spectral energy above the 2/3 dealias cutoff -- the
        under-resolution alarm."""
        p = np.abs(fft(self.omega)) ** 2
        cutoff = self.N // 3
        hi = np.abs(self.k) > cutoff
        return float(np.sum(p[hi]) / (np.sum(p) + 1e-30))

    def diagnostics(self):
        u, u_x = self.velocity()
        wmax = float(np.max(np.abs(self.omega)))
        loc = int(np.argmax(np.abs(self.omega)))
        energy = 0.5 * np.mean(self.omega ** 2) * self.L
        omega_x = self.deriv(self.omega)
        enstrophy = 0.5 * np.mean(omega_x ** 2) * self.L
        return {
            "omega_max": wmax,
            "omega_argmax_x": float(self.x[loc]),
            "ux_max": float(np.max(np.abs(u_x))),
            "energy": float(energy),
            "enstrophy": float(enstrophy),
            "tail_fraction": self.tail_fraction(),
        }


def run_to_time(solver, T, cfl=0.1, dt_max=1e-2, blowup_thresh=1e10,
                tail_thresh=0.05, sample_every=0, on_sample=None):
    """Integrate to time T (or until blowup / under-resolution).

    Returns dict with final time, status, and (if sampling) history list.
    status in {"reached_T", "blowup", "under_resolved"}.
    """
    t = 0.0
    step = 0
    history = []
    status = "reached_T"
    while t < T:
        wmax = solver.omega_max()
        if wmax > blowup_thresh:
            status = "blowup"
            break
        tf = solver.tail_fraction()
        if wmax > 10.0 and tf > tail_thresh:
            status = "under_resolved"
            break
        if sample_every and step % sample_every == 0:
            d = solver.diagnostics()
            d["time"] = t
            history.append(d)
            if on_sample is not None:
                on_sample(t, solver)
        dt = solver.auto_dt(cfl=cfl, dt_max=dt_max)
        dt = min(dt, T - t)
        solver.step_rk4(dt)
        t += dt
        step += 1
    return {"t_final": t, "status": status, "n_steps": step, "history": history,
            "omega_max_final": solver.omega_max()}
