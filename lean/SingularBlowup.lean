/-
# SingularBlowup — Track 1 (algebraic backbone of Conjecture 2.4)

Aggregator root. Every Track-1 file must be imported here so `lake build` (and CI)
actually elaborate it — an un-imported file compiles to nothing (vacuous green).

Scientific context: the 1D Hou–Luo model in self-similar (rescaled) coordinates,
eq. (2.4) of Chen–Huang–Li, arXiv:2604.01868:

    Ω_τ + (U + c_l X) Ω_X = c_ω Ω + Θ_X
    Θ_τ + (U + c_l X) Θ_X = (c_l + 2 c_ω) Θ
    U_X = H(Ω),   U(0) = 0           (Hilbert transform on ℝ)

with the *singular* steady profile (Theorem 2.3):

    Ω̄(X) = 1{X>1}/√(X-1),   Θ̄(X) = (π/2)·1{X>1},   c_l = 2,  c_ω = -1,
    Ū(X) = 2√(1-X) - 2 (X<1),  -2 (X>1).

This repo's numerics (validate/verify_singular_steady.py, src/stability.py, the
co-evolving march) certify the profile and give strong evidence for its linear
stability (Conjecture 2.4). Track 1 here formalises the *closed* algebraic
skeleton of the stability proof. The hard analytic estimate (the (★) weighted-
energy coercivity, L1–L4) is Track 2 and is deliberately NOT stated as a `sorry`.
-/
import SingularBlowup.Profile
import SingularBlowup.Energy
import SingularBlowup.PointSpectrum
import SingularBlowup.SlitSchurMargin
import SingularBlowup.LinearStrictGapClosure
import SingularBlowup.NonlinearComoving
import SingularBlowup.NonlinearBootstrap
import SingularBlowup.CoupledHSector
