# Global Nonlinear Pure-Step ($\eta \equiv 0$) Weak-Edge Compatible Manifold via Autonomous $s$-Clock and Exact Conjugation (2026-11-21)

**Context:** Direct extension of [findings_session_2026_11_21_nonlinear_cusp_manifold.md](findings_session_2026_11_21_nonlinear_cusp_manifold.md) and [findings_session_2026_11_21_pure_step_quadratic_closure.md](findings_session_2026_11_21_pure_step_quadratic_closure.md). Shows that the **exact same autonomous $s$-clock + conjugation method** constructs a full nonlinear global weak-edge compatible manifold **even in the original pure-step ($\eta \equiv 0$) ansatz of the manuscript** (with $\Theta \equiv \frac{\pi}{2}\mathbf{1}_{X>1}$, $B(\tau) \equiv A(\tau)^2 = 1$, $a(\tau) \equiv 0$ for all $\tau \ge 0$, and **no $\sqrt{v}$ $\Theta$-cusp**).

**Evidence tiers:** `[D]` algebraic reduction & asymptotics; `[C]` strip-analytic contraction sketch; `[N]` numerical verification (`validate/edge_pure_step_nonlinear_manifold.py`, `PM-1..PM-3 PASS`, plus independent real-space verification to `2.7e-9` in `findings_session_2026_11_21_weak_edge_hostile_review.md`).

> **Hostile-review framing correction (2026-11-22, A5 & C1):** Working in the conjugated space $\widetilde G = \mathcal{N}_g h$ with $h = O(u^{-2})$ **imposes** the inflow boundary condition $G(0,\tau) = 0$ and weak-edge compatibility $a(\widetilde G(s)) \equiv 0$, $c_1(\widetilde G(s)) = 3a(\text{forcing})$ by choice of function space (it holds for any G-class forcing, not as a property specific to $\mathcal{P}$). Uniqueness of weak solutions on $v > 0$ (where $Q(1) = 0$ allows arbitrary edge-amplitude drifts) is not claimed. What is nontrivial is that the resulting conjugated nonlinear fixed point solves the **real-space** pure-step PDE (independently verified to `2.7e-9` in the core and `3.1e-12` at the edge on `[-104, 24]`, unfiltered).

---

## 1. Exact Conjugated Pure-Step Nonlinear Equation (`[D, N]`)

Setting $\eta \equiv 0$ in the exact autonomous $s$-clock system ($w = \lambda(\tau) v$, $ds = \sqrt{\lambda(\tau)}\,d\tau$, $\lambda(\tau) = \exp(-\int_0^\tau(2 - V_G)\,d\tau')$) yields the single-equation autonomous PDE on $s \in [0, s_\infty)$ ($s_\infty = 1 + O(\varepsilon)$):
$$\widetilde G_s - \frac{1}{2} Q\widetilde G = -(Q\widetilde G)\,D\widetilde G + \frac{1}{2}(Q\widetilde G)\,\widetilde G =: \mathcal{P}(\widetilde G), \qquad D := w\partial_w. \tag{PS}$$
Under the exact Mellin conjugation $\widetilde G = \mathcal{N}_g h$ ($\frac{1}{2}Q(\mathcal{N}_g h) = -\mathcal{N}_r(u h)$, $u = w^{-1/2}$, $N(m) = \frac{\cos(\pi m)\Gamma(2m+2)}{\pi}$), (PS) becomes the forward Duhamel integral equation
$$h(u, s) = e^{-u s} h_0(u) + \int_0^s e^{-u(s - s')}\,\bigl(\mathcal{N}_g^{-1}\mathcal{P}(\mathcal{N}_g h(s'))\bigr)(u)\,ds', \tag{PS-D}$$
equivalently $\widetilde G(s) = S_g(s)(\varepsilon G_0) + \int_0^s S_g(s - s')\,\mathcal{P}(\widetilde G(s'))\,ds'$ with $S_g(s) := \mathcal{N}_g e^{-u s}\mathcal{N}_g^{-1}$.

### Automatic Nonlinear Weak-Edge Compatibility $a(\widetilde G(s)) \equiv 0$ (`[D, N]`, `PM-1`)
Take any small strip-analytic flat initial seed $h_0$ decaying faster than $u^{-1}$ as $u \to \infty$ (e.g. Codex's single-seed profile $h_0(u) = u^3 e^{-u}$, which is holomorphic on every strip $|\operatorname{Im}z| < \rho < \pi$ since $\operatorname{Re}(u) = e^{-\operatorname{Re}(z)/2}\cos(\operatorname{Im}(z)/2) > 0$). Then:
1. For any $\widetilde G \in \mathcal{X}_\rho$ with $a(\widetilde G) = 0$, $\mathcal{P}(\widetilde G) = O(w^2)$ as $w \to 0^+$.
2. Because $N(m)$ has a simple zero at $m = \frac{1}{2}$ with $N'(\frac{1}{2}) = -2$, $\mathcal{N}_g^{-1}\mathcal{P}(\widetilde G)$ has a simple pole at $m = \frac{1}{2}$ with residue $+\frac{\pi}{2}a(\mathcal{P}(\widetilde G))$, hence
   $$f_{\mathcal{P}}(u, s) := \bigl(\mathcal{N}_g^{-1}\mathcal{P}(\widetilde G(s))\bigr)(u) \sim -\frac{\pi}{2}\,a(\mathcal{P}(\widetilde G(s)))\,\frac{1}{u} \qquad (u \to \infty).$$
3. Integrating against $e^{-u(s - s')}$ gains an extra factor of $1/u$ as $u \to \infty$:
   $$h(u, s) = e^{-u s} h_0(u) + \int_0^s e^{-u(s - s')} f_{\mathcal{P}}(u, s')\,ds' \sim -\frac{\pi}{2}\,a(\mathcal{P}(\widetilde G(s)))\,\frac{1}{u^2} = O(u^{-2}) \qquad (u \to \infty).$$
4. Consequently, for **all** $s \in [0, s_\infty)$ (and hence all physical time $\tau \in [0, \infty)$):
   - **Exact weak-edge compatibility:** $a(\widetilde G(s)) = -\frac{2}{\pi}\lim_{u\to\infty} u\,h(u, s) \equiv 0$ (no $\Theta$-cusp needed!).
   - **Exact all-orders nonlinear edge slope:** evaluating the $m = 1$ residue ($N(1) = -6/\pi$) gives
     $$c_1(\widetilde G(s)) = +3\,a(\mathcal{P}(\widetilde G(s))) = -\frac{3}{\pi}\,J(\widetilde G(s)),$$
     extending the $O(\varepsilon^2)$ turning-point formula of [findings_session_2026_11_21_pure_step_quadratic_closure.md](findings_session_2026_11_21_pure_step_quadratic_closure.md) to full nonlinearity without turning-point matching or shooting.

### Simpler Strip-Analytic Contraction (`[D]`)
Because $\mathcal{P}(\widetilde G)$ is already in G-class, $S_g(s - s')$ has **0 derivative loss** and $\mathcal{P}(\widetilde G)$ has **1 derivative loss integrated in $s'$**. Thus (PS-D) is a **single-variable** Nishida–Ovsyannikov contraction on $\mathcal{X}_\rho = H^2_\rho(\frac{1}{4}) \cap H^2_\rho(\frac{3}{4})$ (no $W$ variable, no R-class $\to$ G-class shift, and no $m = \frac{1}{2}$ pole in $h(u,s)$ since $h = O(u^{-2})$!).

*(Adversarial self-audit note on the $\sqrt{v}$-cusp manifold [findings_session_2026_11_21_nonlinear_cusp_manifold.md](findings_session_2026_11_21_nonlinear_cusp_manifold.md) §3: on $\operatorname{Re} m = \sigma_2 > \frac{1}{2}$, the kernel $\frac{u}{u+1}e^{-ut}$ has Mellin norm $\sim t^{-2\sigma_2}$ as $t \to 0^+$, so passing the R-class $\to$ G-class term $\mathcal{T}[\mathcal{R}(\widetilde G)]$ into $H^2_\rho(\sigma_2)$ uses the splitting $\mathcal{I}_0(u,s) = \frac{c_0(s)}{u+1} - \frac{e^{-us}c_0(s)}{u+1} + \int_0^s \frac{u e^{-u(s-s')}}{u+1}(c_0(s') - c_0(s))\,ds'$, which is integrable for $\sigma_2 \in (\frac{1}{2}, 1)$ (e.g. $\sigma_2 = \frac{5}{8}$) because $c_0(s)$ is Lipschitz in $s$, and adding the explicit $O(\varepsilon^2 w^{3/2})$ initial smoothing $h_0(u) = \frac{\pi c_0(0)}{u+1}$ cancels the $s \to 0^+$ term identically. In the pure-step case (PS-D), $h = O(u^{-2})$ has no $m = \frac{1}{2}$ pole at all.)*

---

## 2. Comparison of Both Global Nonlinear Weak-Edge Manifolds (`[D, N]`)

| Property | Pure-Step ($\eta \equiv 0$) Manifold (`PM-1..3`) | $\sqrt{v}$ $\Theta$-Cusp Manifold (`NM-1..3`) |
|---|---|---|
| Weak-edge condition across $X=1$ | $B(\tau) \equiv A(\tau)^2 = 1$, $a(\tau) \equiv 0$, $\eta \equiv 0$ | $B(\tau) \equiv A(\tau)^2 = 1$, $a(\tau) + D(\tau) \equiv 0$ |
| Admissible initial data $G_0$ | Small strip-analytic **edge-flat** $G_0 = \mathcal{N}_g h_0$ ($a(G_0) = 0$) | **Every** small strip-analytic $G_0 \in \mathcal{X}_{\rho_0}$ (including edge-analytic $a(G_0) \neq 0$) |
| Nonlinear shooting needed? | **No** (forward conjugated Duhamel) | **No** (forward conjugated Duhamel, unmodified linear cusp $\eta_0$) |
| Leading edge behavior ($s > 0$) | $G(w,s) = c_1(s)w + \dots$, $c_1 = -\frac{3}{\pi}J(G)$ | $G(w,s) = O(w)$, $\eta(w,s) = \frac{-a_0}{1 + \frac{1}{2}a_0 s}\sqrt{w} + O(w^{3/2})$ |
| Global physical-time decay | $\|G(\cdot,\tau)\|_{L^2_\sigma} \le C\varepsilon e^{-2\sigma\tau}$ | $\|G(\cdot,\tau)\|_{L^2_\sigma} \le C\varepsilon e^{-2\sigma\tau}$ |

---

## 3. Verification Summary (`validate/edge_pure_step_nonlinear_manifold.py`, `PM-1..PM-3`)

- **`PM-1`:** $u f_{\mathcal{P}_1}(u) = +0.310798$ and $u^2 h_2(u) = +0.311223$ at $z = -16$ match $-\frac{\pi}{2} a(\mathcal{P}_1) = +0.310798$ (`rel = 2.2e-6`, `1.4e-3`), with $|u h_2(z=-16)| = 1.04\times 10^{-4}$.
- **`PM-2`:** Full nonlinear pure-step ($\eta \equiv 0$) Picard iteration ($\varepsilon = 0.02$, $s \in [0, 1.2]$, $N_z = 4096, N_s = 240$) converges in 17 iterations (`||dG|| = 7.89e-13`), with $s$-clock PDE residual `6.66e-5` and `max |a(s)| = 2.80e-6`.
- **`PM-3`:** Physical-time reconstruction ($s_\infty = 0.998846$, $\tau_{\max} = 5.68$) satisfies the original physical-time pure-step nonlinear PDE (`res_G(*) = 3.33e-6`) and `max |a(tau)| = 1.17e-7` across $\tau \in \{0.2, 0.6, 1.5, 3.0\}$.
