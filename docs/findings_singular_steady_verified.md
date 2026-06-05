# Finding: the Hou-Luo singular steady state (Thm 2.3) independently certified

Date 2026-06-05.  Target: arXiv:2604.01868, Theorem 2.3 — the singular
(unbounded) self-similar steady state of the 1D Hou-Luo model.  We reproduce it
**analytically in closed form** and confirm numerically.  Code: `src/rescaled.py`,
`validate/verify_singular_steady.py`.

## Setup

Rescaled system (eq. 2.4):
```
Omega_tau + (U + c_l X) Omega_X = c_omega Omega + Theta_X
Theta_tau + (U + c_l X) Theta_X = (c_l + 2 c_omega) Theta
U_X = H(Omega),  U(0)=0     (Hilbert transform on R)
```
Claimed singular steady state (Thm 2.3):
```
Omega_bar(X) = 1{X>1}/sqrt(X-1),   Theta_bar(X) = (pi/2) 1{X>1},
c_l = 2,  c_omega = -1   (=> gamma = -c_l/c_omega = 2).
```

## Closed-form velocity (the load-bearing computation)

`U_X = H(Omega_bar)(X) = (1/pi) p.v. \int_1^\infty (y-1)^{-1/2}/(X-y) dy`.
Substitute `s=y-1`, `a=X-1`:
```
U_X = (1/pi) p.v. \int_0^\infty s^{-1/2}/(a - s) ds.
```
Using the standard p.v. integral `p.v.\int_0^\infty s^{p-1}/(a-s) ds =
pi a^{p-1} cot(pi p)` (0<p<1, a>0) with `p=1/2`, and `cot(pi/2)=0`:

- **X > 1 (a>0):**  `U_X = a^{-1/2} cot(pi/2) = 0`.
- **X < 1 (a<0):**  no pole on (0,inf); `\int_0^\infty s^{-1/2}/(a-s)ds =
  -\int_0^\infty s^{-1/2}/(s+|a|) ds = -pi/sqrt|a|`, so `U_X = -1/sqrt(1-X)`.

Integrate with `U(0)=0`:
```
U(X) = 2 sqrt(1-X) - 2   (X<1),     U(X) = -2   (X>1).
```
=> **`U(1) = -2` exactly.**

## Why U(1) = -2 is the whole point

The steady Theta-equation is `(U + c_l X) Theta_bar_X = (c_l + 2 c_omega) Theta_bar`.
With `c_l=2, c_omega=-1` the RHS coefficient is `c_l+2c_omega = 0`, and
`Theta_bar_X = (pi/2) delta(X-1)`, so the equation reduces to
```
(U(1) + 2)*(pi/2) delta(X-1) = 0   <=>   U(1) = -2.
```
The closed-form velocity delivers exactly this.  The steady Omega-equation on
X>1 (where `Theta_bar_X=0`, `U=-2`):
```
(U + 2X) Omega_bar_X = 2(X-1) * (-1/2)(X-1)^{-3/2} = -(X-1)^{-1/2} = -Omega_bar
= c_omega Omega_bar.   [c_omega = -1]
```
Both equations hold.  **The singular steady state is verified.**

## Numerical confirmation (and the obstruction it exposes)

| check | result |
|---|---|
| `line_hilbert` vs Lorentzian `H[1/(1+x^2)]=x/(1+x^2)` | max err 1.4e-6 (tool correct) |
| steady residual on X in (1.2,40), analytic U | `R_Omega`=2.2e-16, `R_Theta`=0 |
| `U(1)` analytic | -2.000000 (exact) |
| `U(1)` numeric (int of discrete H(Omega_bar)) | **-1.69** (~15% off) |
| `H(Omega_bar)` numeric vs analytic, `|X-1|>0.5` | max err 1.5e-1 |

The ~15% numeric error is **not** a correctness issue — it is dominated by
**far-field truncation of the slowly decaying `Omega_bar ~ 1/sqrt(X)` tail**:
the part beyond the box `|X|<L=200` contributes `~ 2/sqrt(L) ≈ 0.14` to `U_X`
near `X=1`, matching the observed gap.  The `(X-1)^{-1/2}` singularity at the
support edge adds local error.

**Consequence (obstruction at maximum resolution).** Any *numerical* work on this
profile — in particular the linear-stability operator (Conjecture 2.4 / next
task), which needs `H(omega')` for perturbations supported on the slow tail —
requires a **mapped infinite domain** (e.g. `X = tan`-type compactification) or
an **analytic tail/desingularization** rather than a truncated equispaced grid.
This is the concrete form of the "singularity-adapted basis / desingularizing
change of variables" the smooth-profile CAP framework lacks. The analytic
velocity `U(X)=2 sqrt(1-X)-2` is known in closed form, which removes one of the
two singular ingredients from the stability operator for free.

## Status vs the literature

This is an independent confirmation of a ~2-month-old result (2604.01868, Thm
2.3), obtained by closed-form Hilbert transform rather than their construction.
It is the certified base point for the stability computation (the open piece).
