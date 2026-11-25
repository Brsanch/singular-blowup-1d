"""Exact singular self-similar profiles and linear stability operators for gCLM (arXiv:2603.25104)."""

from .profiles import (
    mu_of_a,
    scaling_exponents,
    omega_bar,
    omega_bar_X,
    Ux_bar,
    U_bar,
    advection_bar,
    steady_residual,
    interior_stagnation_point,
    hilbert_omega_bar_quad,
    mellin_symbol_c,
    conjugation_N_am1,
)
from .operator import (
    build_operator_z,
    build_operator_zeta,
    nonlinear_rhs_z,
    spectrum,
)

__all__ = [
    "mu_of_a",
    "scaling_exponents",
    "omega_bar",
    "omega_bar_X",
    "Ux_bar",
    "U_bar",
    "advection_bar",
    "steady_residual",
    "interior_stagnation_point",
    "hilbert_omega_bar_quad",
    "mellin_symbol_c",
    "conjugation_N_am1",
    "build_operator_z",
    "build_operator_zeta",
    "nonlinear_rhs_z",
    "spectrum",
]
