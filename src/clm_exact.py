#!/usr/bin/env python3
"""
Exact closed-form solution of the Constantin-Lax-Majda model (a=0, nu=0):

    omega_t = omega * H[omega].

Derivation (self-contained, so the gate certifies itself):
Let a = H[omega], b = omega.  The Tricomi identity H(2 omega H omega) =
(H omega)^2 - omega^2 gives, after applying H to the PDE,
    a_t = (1/2)(a^2 - b^2),   b_t = a b.
Hence  w = a + i b = H[omega] + i*omega  obeys the pointwise Riccati ODE
    w_t = (1/2) w^2     =>     w(x,t) = w0(x) / (1 - (t/2) w0(x)),
with w0 = H[omega0] + i*omega0.  Taking the imaginary part:

    omega(x,t) = 4 * omega0 / ( (2 - t*H[omega0])^2 + t^2 * omega0^2 ).      (*)

Blowup: the denominator vanishes where omega0(x*)=0 and H[omega0](x*)>0, at
    t*(x*) = 2 / H[omega0](x*),
so the first singular time is  t* = 2 / max{ H[omega0](x) : omega0(x)=0 }.

For omega0 = sin(x) on [0,2*pi):  H[sin] = -cos, zeros at x=0,pi with
H[omega0] = -1, +1 -> first blowup at x*=pi, t*=2.  omega itself stays 0 at the
center and blows up like 1/xi just off it (max|omega| ~ 1/(2-t)).
"""

import numpy as np


def clm_exact(solver, omega0, t):
    """Exact CLM field at time t, using the solver's Hilbert convention.

    omega0 : initial vorticity on solver.x
    returns omega(x, t) on the same grid via formula (*).
    """
    Hw0 = solver.hilbert(omega0)
    denom = (2.0 - t * Hw0) ** 2 + (t ** 2) * (omega0 ** 2)
    return 4.0 * omega0 / denom


def clm_blowup_time(solver, omega0):
    """First singular time t* = 2 / max{ H[omega0] where omega0 = 0 }.

    Returns (t_star, x_star).  Found by locating sign changes of omega0 and
    interpolating H[omega0] at the zero crossing (so it is grid-robust).
    """
    Hw0 = solver.hilbert(omega0)
    x = solver.x
    N = solver.N
    best = (np.inf, None)
    for i in range(N):
        a, b = omega0[i], omega0[(i + 1) % N]
        if a == 0.0:
            Hz = Hw0[i]
            xz = x[i]
        elif a * b < 0.0:  # sign change between i and i+1
            frac = a / (a - b)  # linear zero location in [0,1]
            Hz = Hw0[i] + frac * (Hw0[(i + 1) % N] - Hw0[i])
            xz = x[i] + frac * solver.dx
        else:
            continue
        if Hz > 0.0:
            ts = 2.0 / Hz
            if ts < best[0]:
                best = (ts, xz)
    return best
