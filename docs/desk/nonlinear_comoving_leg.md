# Direction 2: Nonlinear Co-Moving Leg (Exact Gauge-(a) Collapse & Dual-Strip Architecture)

**Date:** 2026-11-20
**Status:** Exact nonlinear PDE collapse, dual-strip $(\sigma_+ = \frac{1}{4}, \sigma_- = \frac{3}{4})$ symbol/resolvent closure, and Hilbert edge commutator **CLOSED & VERIFIED** (`validate/edge_nonlinear_comoving_closure.py` `NL-1..NL-5 PASS`; Lean 4 `lean/SingularBlowup/NonlinearComoving.lean` `0 sorry`, `0 axioms`, vacuity lint `0 findings`). Remaining continuous PDE step: full $H^2_{\mathrm{edge}}([1,2])$ quasilinear Duhamel/energy bootstrap closing the $Y$-coordinate edge regularity.

---

## 1. Exact Gauge-(a) Co-Moving Collapse (Zero Approximation)

Starting from the dynamic-rescaling 1D Hou–Luo system (Chen–Huang–Li, arXiv:2604.01868 eq. (2.4)) on $X > 1$ with $z = \ln(X - 1)$:
$$\Omega(X, \tau) = (X - 1)^{-1/2}\bigl(1 + G(z, \tau)\bigr)\mathbf{1}_{X>1}, \qquad \Theta(X, \tau) = B(\tau)\left(\frac{\pi}{2} + h(z, \tau)\right)\mathbf{1}_{X>1},$$
impose **Gauge (a)** ($A(\tau) \equiv 1$, edge pinned at $X = 1$):
$$c_l(\tau) = -U(1, \tau) = 2 - V_G(1), \qquad c_\omega(\tau) = -1 + \frac{V_G(1) - H[g](1)}{2}, \qquad \frac{\dot{B}}{B} = -H[g](1),$$
where:
$$V_G(1) = -\frac{1}{\pi}\int_{-\infty}^\infty e^{z/2}\ln(1 + e^{-z})\,G(z)\,dz, \qquad H[g](1) = -\frac{1}{\pi}\widehat{G}\!\left(\frac{1}{2}\right) = -\frac{1}{\pi}\int_{-\infty}^\infty e^{-z/2} G(z)\,dz.$$

Three exact cancellations occur simultaneously (Lean: `gauge_a_advection_velocity_collapse`, `gauge_a_exact_omega_collapse`, `gauge_a_exact_theta_balance`; Gate `NL-1`: `err_G = 4.16e-12`, `err_h = 2.22e-16`):
1. **Zero singular $e^{-z}$ dilation drift:** $U(Y) + c_l Y = (Y - 1)\bigl(2 - V_G(1) + Q_{\mathrm{op}} G(z)\bigr)$.
2. **Zero inhomogeneous amplitude forcing in $G$:** Setting $\widetilde{Q}_{\mathrm{op}} G := Q_{\mathrm{op}} G - H[g](1)$, all constant and linear $(V_G(1), H[g](1))$ terms cancel identically.
3. **Zero inhomogeneous $\Theta$-forcing in $h$:** Since $c_l + 2 c_\omega = -H[g](1) = \dot{B}/B$, $h$ solves pure nonlinear transport.

The **exact, zero-approximation** nonlinear co-moving PDE is:
$$\begin{aligned}
\partial_\tau G &= \underbrace{-2\partial_z G + \tfrac{1}{2}\widetilde{Q}_{\mathrm{op}} G}_{\widetilde{M}_{GG} G} + B(\tau) e^{-z/2}\partial_z h + \bigl(V_G(1) - Q_{\mathrm{op}} G\bigr)\partial_z G + \tfrac{1}{2}\bigl(\widetilde{Q}_{\mathrm{op}} G\bigr) G, \\
\partial_\tau h &= -2\partial_z h + \bigl(V_G(1) - Q_{\mathrm{op}} G\bigr)\partial_z h, \\
\dot{B} &= -B\,H[g](1).
\end{aligned}$$

---

## 2. Dual-Strip $(\sigma_+ = \frac{1}{4}, \sigma_- = \frac{3}{4})$ Structure & Exact Rational Symbol Bounds

On the two-weight intersection space $L^2_{3/4, 1/4}(\mathbb{R}) = L^2_{3/4}(\mathbb{R}) \cap L^2_{1/4}(\mathbb{R})$:
- **Exact unit-tangent identity:** $|\tan(\pi(\frac{1}{4} + i\tau))| = |\tan(\pi(\frac{3}{4} + i\tau))| \equiv 1$ for all $\tau \in \mathbb{R}$, yielding exact rational operator norms for $T_c$ ($c(s) = -\frac{\tan(\pi s)}{2s+1}$) and $T_c \partial_z$ (Lean: `symbol_norm_quarter_sq_le`..`lower_order_smoothing_const`; Gate `NL-2`):
  $$\|T_c\|_{L^2_{1/4}} = \frac{2}{3}, \qquad \|T_c\|_{L^2_{3/4}} = \frac{2}{5}, \qquad \|T_c \partial_z\|_{L^2_{1/4}} = \|T_c \partial_z\|_{L^2_{3/4}} = \frac{1}{2},$$
  $$|H[g](1)| \le \frac{1}{\sqrt{\pi}}\|G\|_{L^2_{3/4, 1/4}}, \qquad |V_G(1)| \le 0.5104\,\|G\|_{L^2_{1/4}},$$
  and the **lower-order derivative smoothing bound** (zero top-order derivative on the nonlocal term!):
  $$\left\|\frac{1}{2}\partial_z(\widetilde{Q}_{\mathrm{op}} G)\right\|_{L^2_{1/4}} \le \frac{7}{10}\,\|G\|_{L^2_{3/4}}.$$
- **Dual-strip `W1` symmetry:** $L(\sigma, \delta) = 3\pi + \frac{2}{\delta} + \frac{1}{\sigma} + \frac{1}{1-\sigma}$ satisfies $L(\frac{3}{4}, \delta) = L(\frac{1}{4}, \delta)$ (`dual_strip_digamma_symmetry`; Gate `NL-3`: `W1` Route F holds at $\sigma_- = 3/4, \delta = 1/20$, $\operatorname{Re}\lambda = -1.45$).
- **Automatic $s = 1/2$ vanishing on $D(M_{GG})$:** In the W1 slit-plane resolvent $\widehat{G}(s;\lambda) = X(s;\lambda)[D(e^{4\pi i s}) - D(1)]$, evaluating at $s = 1/2$ gives $e^{4\pi i(1/2)} = 1 \implies \widehat{G}(1/2;\lambda) = 0 \implies H[g](1) = 0$ identically on the linear resolvent image, matching $\widetilde{M}_{GG} = M_{GG}$ on $D(M_{GG})$.

---

## 3. Edge Hilbert Commutator & Remaining Nonlinear PDE Step

- **Why fractional edge weights lose $\frac{1}{2}$ while smooth $Y$-perturbations preserve exponent $1$:**
  Because $c(s) = -\frac{\tan(\pi s)}{2s+1}$ vanishes at $s = 1$, for smooth $Y$-perturbations $F(Y) = (Y-1) F_1(Y)$ ($\alpha = 1$), the exact fraction identity $\frac{t-1}{Y-t} = \frac{Y-1}{Y-t} - 1$ (`hilbert_edge_commutator_fraction`) gives
  $$H\bigl[(t-1)^{-1/2} F\bigr](Y) - H[g](1) = (Y - 1)\, H\bigl[(t-1)^{-1/2} F_1\bigr](Y),$$
  preserving $\widetilde{Q}_{\mathrm{op}} F(Y) = O(Y - 1)$ with **zero exponent loss** (Gate `NL-4`: `slope(alpha=1.00) = 1.0001` vs `slope(alpha=0.75) = 0.2505`).
- **Self-consistent nonlinear verification (Gate `NL-5`):** Under the exact Gauge-(a) flow (no ad-hoc projection), edge-smooth perturbations relax monotonically in every norm (`||G||_{3/4,1/4}: 0.312 -> 0.0345`, `||G_z||_{3/4,1/4}: 0.270 -> 0.0256`, `||Q_op G||_inf: 0.141 -> 0.0035`, `||F_1||_{L2[1,2]}: 0.164 -> 2.02e-5`, `c_l -> 2.0029`, `c_om -> -1.00001`, `B -> 1.1682`).
- **Named remaining continuous PDE obligation for complete Theorem 2.4:**
  Close the $H^2([1,2])$ energy estimate for $F_1(Y) = F(Y)/(Y-1)$ (coupling the physical edge transport damping $-(k + \frac{1}{2}) a_1$ with the bounded $L^2([1,2])$ Hilbert commutator $H[(t-1)^{1/2} F_2]$) and stitch it to the far-field $H^2_{1/4}([0,\infty))$ quasilinear energy estimate (`nonlinear_transport_damping_margin`) + `W1–W3` linear semigroup Duhamel contraction.
