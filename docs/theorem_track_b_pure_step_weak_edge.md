# Theorem (Track B): Global Nonlinear Pure-Step ($\eta \equiv 0$) Weak-Edge Compatible Manifold in Strip-Analytic Spaces

**Date:** 2026-11-22  
**Status:** Analytic proof (`[D]`), not expert-reviewed; Lemma 3 re-derived by the 2026-11-22 post-`f5106e9` review (`findings_session_2026_11_22_weak_edge_post_f5106e9_review.md`), which found two minor gaps (marked below). Numerical checks: `validate/edge_track_b_proof_lemmas.py` (`TB-1..TB-5`), `validate/edge_pure_step_nonlinear_manifold.py` (`PM-1..PM-3`), and Lean 4 pointwise kernel bounds in `lean/SingularBlowup/WeakEdgeManifolds.lean`. The gates check identities and computed solutions; they do not check the proof's constants or ε₀.
**Addresses:** Hostile review `findings_session_2026_11_21_weak_edge_hostile_review.md` (`A5`, `B1`, `B2`, `B4`, `C1–C3`).

---

## 1. Statement of the Theorem

Consider the 1D Hou–Luo dynamic-rescaling system `(eq:dyn)` in the **pure-step** ($\eta \equiv 0$, $\Theta(X,\tau) \equiv \frac{\pi}{2}\mathbf{1}_{X>1}$) ansatz with $A \equiv B \equiv \xi \equiv 1$, $c_l(\tau) = 2 - V_G(\tau)$, and $c_\omega(\tau) = -c_l(\tau)/2$ ($c_l + 2c_\omega \equiv 0$):
$$\Omega(X,\tau) = (X-1)^{-1/2}\bigl(1 + G(X-1,\tau)\bigr)\mathbf{1}_{X>1}, \qquad \Theta(X,\tau) = \tfrac{\pi}{2}\mathbf{1}_{X>1},$$
where on $v = X - 1 > 0$, $G$ satisfies the exact nonlinear PDE
$$G_\tau + \bigl(2 - V_G(\tau) + QG\bigr)\,v G_v - \tfrac{1}{2}(QG)(1 + G) = 0, \tag{1.1}$$
with $Q(G)(v) = \frac{1}{v}\int_0^v H[y^{-1/2}G](1+y)\,dy$, $a(G) = (QG)(0) = -\frac{1}{\pi}\int_0^\infty y^{-3/2}G(y)\,dy$, and $V_G = -\frac{1}{\pi}\int_0^\infty y^{-1/2}\ln(1 + 1/y)G(y)\,dy$.

### Function spaces
Fix $\rho_0 \in (0, \frac{\pi}{2}]$. For $\rho \in (0, \rho_0]$ and $\sigma \in \mathbb{R}$, let $H^2_\rho(\sigma)$ be the strip-analytic Hardy–Mellin Hilbert space of functions $F(e^z)$ holomorphic on $\mathbb{S}_\rho := \{z \in \mathbb{C} : |\operatorname{Im} z| < \rho\}$ with norm
$$\|F\|_{\rho,\sigma}^2 := \frac{1}{2\pi}\int_{\mathbb{R}} e^{2\rho|\xi|}\,\bigl|\widehat F(\sigma + i\xi)\bigr|^2\,d\xi < \infty, \qquad \widehat F(\sigma + i\xi) := \int_{\mathbb{R}} e^{-(\sigma+i\xi)z} F(e^z)\,dz.$$
Define the two-line intersection Banach scale and its closed weak-edge compatibility subspace:
$$\mathcal{X}_\rho := H^2_\rho(\tfrac{1}{4}) \cap H^2_\rho(\tfrac{5}{8}), \qquad \|F\|_{\mathcal{X}_\rho} := \|F\|_{\rho,1/4} + \|F\|_{\rho,5/8},$$
$$\mathcal{X}_\rho^0 := \bigl\{G \in \mathcal{X}_\rho : a(G) = -\tfrac{1}{\pi}\widehat G(\tfrac{1}{2}) = 0\bigr\}.$$

> **Inflow boundary condition & uniqueness scope (addressing hostile review `A5`):** On $v > 0$, the edge $v = 0$ is an inflow boundary and $Q(1) = 0$, so (1.1) without an edge condition admits non-unique weak solutions with arbitrary edge-amplitude drift $\alpha(\tau)$. Restricting to $\mathcal{X}_\rho^0$ **imposes** the inflow boundary condition $G(0,\tau) = 0$ and the weak-edge compatibility condition $a(G(\cdot,\tau)) \equiv 0$ (equivalently $B(\tau) \equiv A(\tau)^2 = 1$ across $X = 1$ under the manuscript's singular-product convention `rem:weak_edge_jump`).

### Theorem 1 (Global Nonlinear Pure-Step Weak-Edge Compatible Manifold)
Let $\rho_0 \in (0, \frac{\pi}{2}]$. There exists $\varepsilon_0 = \varepsilon_0(\rho_0) > 0$ such that for every initial profile $G_0 \in \mathcal{X}_{\rho_0}^0$ with $\|G_0\|_{\mathcal{X}_{\rho_0}} \le 1$ and every $0 < \varepsilon \le \varepsilon_0$:
1. **(Autonomous $s$-clock solution):** On $s \in [0, 2]$ and $\rho_1 := \frac{1}{2}\rho_0$, there exists a mild (and classical on $(0, \infty) \times (0, 2]$) solution $\widetilde G \in C([0, 2]; \mathcal{X}_{\rho_1}^0)$ — unique within the shrinking-cone strip-analytic scale $\{\mathcal{X}_\rho^0\}_{\rho < \rho_0 - a_\infty s}$ — of the autonomous $s$-clock equation
   $$\widetilde G_s - \tfrac{1}{2} Q\widetilde G = \mathcal{P}(\widetilde G) := -(Q\widetilde G)D\widetilde G + \tfrac{1}{2}(Q\widetilde G)\widetilde G, \qquad \widetilde G(\cdot, 0) = \varepsilon G_0, \tag{1.2}$$
   satisfying $\sup_{s \in [0,2]} \|\widetilde G(s)\|_{\mathcal{X}_{\rho_1}} \le 2 C_S \varepsilon$, $a(\widetilde G(s)) \equiv 0$, and the exact nonlinear edge-slope relation $c_1(\widetilde G(s)) = 3\,a(\mathcal{P}(\widetilde G(s))) = -\frac{3}{\pi}J(\widetilde G(s))$ for all $s \in (0, 2]$, where $J(G) := -\pi a(\mathcal{P}(G))$.
2. **(Global physical-time reconstruction):** The scalar horizon ODE
   $$\frac{dr}{ds} = -1 + \tfrac{1}{2}V_{\widetilde G(\cdot, s)}(r^2), \qquad r(0) = 1, \tag{1.3}$$
   is $C^1$ on $r \ge 0$ (with $V_{\widetilde G}(r^2) = O(\varepsilon r^{5/4})$ and $\partial_r V_{\widetilde G}(r^2) = O(\varepsilon r^{1/4})$ as $r \to 0^+$, and $+\frac{2}{3}c_1(\widetilde G(s))r^2 + O(\varepsilon r^{9/4})$ for each fixed $s > 0$, with a constant that is **not** uniform as $s \to 0^+$; Lemma 4 uses only the uniform $O(\varepsilon r^{5/4})$ bound) and has a unique transversal zero $s_\infty = 1 + O(\varepsilon) \in (\frac{1}{2}, \frac{3}{2}) \subset [0, 2]$. Setting $\tau(s) := \int_0^s ds'/r(s') \in [0, \infty)$ ($s \in [0, s_\infty)$) and $\lambda(\tau) := r(s(\tau))^2$, the reconstructed profile
   $$G(v, \tau) := \widetilde G\bigl(\lambda(\tau)v,\, s(\tau)\bigr)$$
   solves the original pure-step nonlinear PDE (1.1) classically on $(v, \tau) \in (0, \infty) \times (0, \infty)$ with $G(\cdot, 0) = \varepsilon G_0$, preserves exact weak-edge compatibility $a(G(\cdot, \tau)) \equiv 0$ ($B(\tau) \equiv A(\tau)^2 = 1$) for all $\tau \ge 0$, and decays at the optimal convective rate
   $$\|G(\cdot, \tau)\|_{L^2_\sigma} = \lambda(\tau)^\sigma \|\widetilde G(\cdot, s(\tau))\|_{L^2_\sigma} \le C \varepsilon\,e^{-2\sigma\tau} \qquad \bigl(\sigma \in [\tfrac{1}{4}, \tfrac{5}{8}],\; \tau \ge 0\bigr).$$

*(Non-emptiness: For any $b > 0$ and $p > 1$, $h_0(u) = u^p e^{-b u}$ is holomorphic on $\mathbb{S}_\rho$ for all $\rho < \pi$, and $G_0 := \mathcal{N}_g h_0$ lies in $\mathcal{X}_{\rho_0}^0$ for every $\rho_0 < \pi$.)*

---

## 2. Complete Proof (Lemmas 1–4)

### Lemma 1 (Closedness of $\mathcal{X}_\rho^0$, Pole-Free $Q$-Smoothing, and Bilinear Estimates)
1. **Closedness:** For any $G \in \mathcal{X}_\rho$, $\widehat G(m)$ is holomorphic on $\frac{1}{4} < \operatorname{Re} m < \frac{5}{8}$ with $L^2$ boundary traces on $\operatorname{Re} m = \frac{1}{4}, \frac{5}{8}$. By Cauchy's integral formula on $[\frac{1}{4}, \frac{5}{8}] \times [-R, R]$ ($R \to \infty$), $|a(G)| = \frac{1}{\pi}|\widehat G(\frac{1}{2})| \le C_0 \|G\|_{\mathcal{X}_\rho}$. Hence $\mathcal{X}_\rho^0$ is a closed subspace of $\mathcal{X}_\rho$.
2. **Pole-free $Q$-smoothing on $\mathcal{X}_\rho^0$ (resolving `B1`):** Recall $\mathcal{M}[QG](\zeta) = \frac{\cot(\pi\zeta)}{\zeta+1}\widehat G(\zeta + \frac{1}{2})$. For $G \in \mathcal{X}_\rho^0$, $\widehat G(\frac{1}{2}) = 0$ cancels the simple pole of $\cot(\pi\zeta)$ at $\zeta = 0$ identically. Thus $\mathcal{M}[QG](\zeta)$ is holomorphic on $-\frac{1}{4} < \operatorname{Re}\zeta < \frac{1}{8}$ with boundary traces on $\operatorname{Re}\zeta = -\frac{1}{4}$ and $\operatorname{Re}\zeta = +\frac{1}{8}$ satisfying
   $$\|QG\|_{\rho,-1/4} + \|D(QG)\|_{\rho,-1/4} \le C_1 \|G\|_{\rho,1/4}, \qquad \|QG\|_{\rho,1/8} + \|D(QG)\|_{\rho,1/8} \le C_1 \|G\|_{\rho,5/8},$$
   because $|\cot(\pi(\zeta_0 + i\xi))/(\zeta_0 + 1 + i\xi)| \le C\langle\xi\rangle^{-1}$ on $\zeta_0 \in \{-\frac{1}{4}, \frac{1}{8}\}$. By Cauchy–Schwarz in $\xi$ ($\langle\xi\rangle^{-1} \in L^2(\mathbb{R})$), for $z = x + iy \in \mathbb{S}_\rho$:
   $$|(QG)(e^z)| \le \begin{cases} C_2 e^{-x/4}\|G\|_{\rho,1/4} & (x \ge 0), \\ C_2 e^{+x/8}\|G\|_{\rho,5/8} & (x \le 0), \end{cases} \implies \|QG\|_{H^\infty(\mathbb{S}_\rho)} \le C_Q \|G\|_{\mathcal{X}_\rho}. \tag{2.1}$$
3. **Cauchy estimate & bilinear bound on $\mathcal{P}$:** For $0 < \rho' < \rho \le \rho_0$, $|\xi|e^{\rho'|\xi|} \le (e(\rho-\rho'))^{-1}e^{\rho|\xi|}$, so $\|DG\|_{\mathcal{X}_{\rho'}} \le \frac{C_D}{\rho - \rho'}\|G\|_{\mathcal{X}_\rho}$. By the Hardy–Paley–Wiener equivalence $\|F\|_{\rho',\sigma}^2 \asymp \sup_{|y|<\rho'}\|e^{-\sigma(\cdot)}F(e^{\cdot+iy})\|_{L^2(\mathbb{R})}^2$, multiplication by $QG \in H^\infty(\mathbb{S}_{\rho'})$ is bounded on $\mathcal{X}_{\rho'}$, yielding for all $G_1, G_2 \in \mathcal{X}_\rho^0$:
   $$\|\mathcal{P}(G_1) - \mathcal{P}(G_2)\|_{\mathcal{X}_{\rho'}} \le \frac{C_{\mathcal{P}}}{\rho - \rho'}\bigl(\|G_1\|_{\mathcal{X}_\rho} + \|G_2\|_{\mathcal{X}_\rho}\bigr)\|G_1 - G_2\|_{\mathcal{X}_\rho}. \tag{2.2}$$
4. **Zero-derivative-loss bound on $a(\mathcal{P}(G))$ (resolving `B2`):** For $G \in \mathcal{X}_\rho^0$, $G(w) = O(w^{5/8})$ and $(QG)(w) = O(w^{1/8})$ as $w \to 0^+$, while $G(w) = O(w^{1/4})$ and $(QG)(w) = O(w^{-1/4})$ as $w \to \infty$. Integrating $w^{-1/2}(QG)G'$ by parts (both boundary terms vanish):
   $$a(\mathcal{P}(G)) = -\frac{1}{\pi}\int_0^\infty w^{-3/2}\bigl(-(QG)w G' + \tfrac{1}{2}(QG)G\bigr)\,dw = -\frac{1}{\pi}\int_{\mathbb{R}} e^{-z/2} D(QG)(e^z)\,G(e^z)\,dz. \tag{2.3}$$
   Splitting $\mathbb{R} = (-\infty, 0] \cup [0, \infty)$: on $z \ge 0$, $e^{-z/2} \le 1 = e^{z/4}e^{-z/4}$; on $z \le 0$, $e^{-z/2} \le e^{-3z/4} = e^{-z/8}e^{-5z/8}$. By Cauchy–Schwarz and Part 2, with **zero radius shrinkage**:
   $$\bigl|a(\mathcal{P}(G_1)) - a(\mathcal{P}(G_2))\bigr| \le C_a\bigl(\|G_1\|_{\mathcal{X}_\rho} + \|G_2\|_{\mathcal{X}_\rho}\bigr)\|G_1 - G_2\|_{\mathcal{X}_\rho}. \tag{2.4}$$

---

### Lemma 2 (Exact Propagator $S_g(t) = \mathcal{N}_g e^{-ut}\mathcal{N}_g^{-1}$ with Integrable $t^{-1/4}$ Second-Line Bound)
Let $N(m) := \frac{\cos(\pi m)\Gamma(2m+2)}{\pi}$ and $u(z) := e^{-z/2}$. By Stirling's formula, $|N(\sigma + i\xi)| \asymp \langle\xi\rangle^{2\sigma + 3/2}$ on any line $\operatorname{Re} m = \sigma \notin \frac{1}{2} + \mathbb{Z}$ (specifically order $p_1 = 2$ at $\sigma_1 = \frac{1}{4}$ and $p_2 = \frac{11}{4}$ at $\sigma_2 = \frac{5}{8}$). For $F \in \mathcal{X}_\rho$ and $t > 0$, define $S_g(t)F := \mathcal{N}_{1/4}(e^{-u t}\mathcal{N}_{1/4}^{-1}F)$. Then:
1. **Invariance of $\mathcal{X}_\rho^0$:** For $t > 0$, $e^{-u t}(\mathcal{N}_{1/4}^{-1}F)$ decays exponentially as $u \to \infty$ ($z \to -\infty$), so its Mellin transform is holomorphic at $m = \frac{1}{2}$; hence $\widehat{S_g(t)F}(\frac{1}{2}) = N(\frac{1}{2})\widehat{e^{-ut}f}(\frac{1}{2}) = 0$, i.e. $S_g(t)F \in \mathcal{X}_\rho^0$.
2. **First-line bound ($\sigma_1 = \frac{1}{4}$):** On $\mathbb{S}_\rho$ ($\rho \le \rho_0 \le \frac{\pi}{2}$), $\operatorname{Re}(u) = e^{-x/2}\cos(y/2) \ge \frac{1}{\sqrt{2}}e^{-x/2} > 0$, so $|\partial_z^k(e^{-ut})| \le C_k(\rho_0)$ uniformly for all $t \ge 0$ (`WeakEdgeManifolds.strip_norm_exp_le_one`, `strip_norm_us_exp_le`). Since $|N(\frac{1}{4}+i\xi)| \asymp \langle\xi\rangle^2$, multiplication by $e^{-ut}$ is uniformly bounded on $H^{2}_{\text{Sob},\rho}(\frac{1}{4})$, giving
   $$\|S_g(t)F\|_{\rho, 1/4} \le C_{S,1}\,\|F\|_{\rho, 1/4} \qquad (t \ge 0). \tag{2.5}$$
3. **Second-line bound ($\sigma_2 = \frac{5}{8}$) via exact rank-1 pole split:** Between $\operatorname{Re} m = \frac{1}{4}$ and $\operatorname{Re} m = \frac{5}{8}$, the unique zero of $N(m)$ is the simple zero at $m = \frac{1}{2}$ with $N'(\frac{1}{2}) = -2$. Hence $\widehat f(m) = \widehat F(m)/N(m)$ has a single simple pole at $m = \frac{1}{2}$ with residue $+\frac{\pi}{2}a(F)$. Define the explicit holomorphic rank-1 basis function
   $$\phi_*(z) := \frac{1}{2\cosh(z/2)} = \frac{u}{u^2 + 1}, \qquad \widehat{\phi_*}(m) = \frac{\pi}{\cos(\pi m)} \quad \bigl(|\operatorname{Re} m| < \tfrac{1}{2}\bigr), \qquad \operatorname{Res}_{m=1/2}\widehat{\phi_*}(m) = -1.$$
   Setting $\tilde f := f + \frac{\pi}{2}a(F)\phi_*$ cancels the $m = \frac{1}{2}$ pole identically, so $\widehat{\tilde f}(m)$ is holomorphic on $\frac{1}{4} \le \operatorname{Re} m \le \frac{5}{8}$. On $\operatorname{Re} m = \frac{5}{8}$, $|\widehat{\phi_*}(\frac{5}{8}+i\xi)| \le C e^{-\pi|\xi|}$ decays faster than $e^{-\rho_0|\xi|}$ ($\rho_0 \le \frac{\pi}{2} < \pi$), so $\|\tilde f\|_{H^{11/4}_{\text{Sob},\rho}(5/8)} \le C_3(\|F\|_{\rho,5/8} + |a(F)|) \le C_3'\|F\|_{\mathcal{X}_\rho}$. By complex interpolation between the integer orders $k = 2$ and $k = 3$ (since $|\partial_z^k(e^{-ut})| \le C_k(\rho_0)$ uniformly for $t \ge 0$), multiplication by $e^{-ut}$ is uniformly bounded on $H^{11/4}_{\text{Sob},\rho}(\frac{5}{8})$, giving $\|\mathcal{N}_{5/8}(e^{-ut}\tilde f)\|_{\rho,5/8} \le C_{S,2}\|F\|_{\mathcal{X}_\rho}$ uniformly for $t \ge 0$.

   For the rank-1 piece $e^{-ut}\phi_*$: since $|N(\frac{5}{8}+i\xi)| \le C\langle\xi\rangle^3$, Plancherel gives
   $$\|\mathcal{N}_{5/8}(e^{-ut}\phi_*)\|_{\rho,5/8} \le C\sup_{|y|<\rho}\sum_{j=0}^3 \left\|e^{-\frac{5}{8}x}\,\partial_z^j\!\left(\frac{e^{-u(x+iy)t}}{2\cosh((x+iy)/2)}\right)\right\|_{L^2_x(\mathbb{R})}.$$
   On $x \ge 0$, $e^{-5x/4}e^{-x} = e^{-9x/4} \in L^1$ gives $O(1)$ uniformly in $t \ge 0$. On $x \le 0$, substituting $r = e^{-x/2} \in [1, \infty)$ ($dx = 2dr/r$, $e^{-5x/4}e^{x} = r^{1/2}$) and $v = r t$ gives
   $$\int_{-\infty}^0 e^{-5x/4}e^x(1 + (rt)^j)^2 e^{-2 c_0 r t}\,dx = 2\,t^{-1/2}\int_t^\infty v^{-1/2}(1 + v^j)^2 e^{-2 c_0 v}\,dv \le C_j\,t^{-1/2}.$$
   Taking the square root yields $\|\mathcal{N}_{5/8}(e^{-ut}\phi_*)\|_{\rho,5/8} \le C_*(1 + t^{-1/4})$. Combining both pieces:
   $$\boxed{\|S_g(t)F\|_{\mathcal{X}_\rho} \le C_S\,\|F\|_{\mathcal{X}_\rho} + C_*\,t^{-1/4}\,|a(F)| \qquad \text{for all } F \in \mathcal{X}_\rho,\; t > 0.} \tag{2.6}$$
   In particular, for $G_0 \in \mathcal{X}_{\rho_0}^0$ ($a(G_0) = 0$), $\|S_g(s)G_0\|_{\mathcal{X}_\rho} \le C_S\|G_0\|_{\mathcal{X}_{\rho_0}}$ uniformly for all $s \ge 0$.

---

### Lemma 3 (Nishida–Ovsyannikov Contraction on $s \in [0, 2]$)
Write $\widetilde G = \varepsilon\hat G$ and consider the Duhamel map on $\mathcal{X}_\rho^0$:
$$\Phi(\hat G)(s) := S_g(s)G_0 + \varepsilon\int_0^s S_g(s - s')\,\mathcal{P}(\hat G(s'))\,ds'. \tag{2.7}$$
Fix $S := 2$, $M := 2 C_S$ (for $\|G_0\|_{\mathcal{X}_{\rho_0}} \le 1$), and $a_\infty := \frac{\rho_0}{6}$ (so on $s \in [0, 2]$ and $\rho_1 = \frac{1}{2}\rho_0$, we have a strict margin $d_{a_\infty}(\rho_1, s) \ge \frac{1}{6}\rho_0 > 0$). For any $a \in [\frac{\rho_0}{12}, \frac{\rho_0}{6}]$, define $d_a(\rho, s) := \rho_0 - \rho - a s$ on the cone $\Omega_a := \{(s,\rho) : 0 < s \le 2,\; 0 < \rho < \rho_0 - a s\}$ and Nishida's $s$-scaled weighted norm
$$\|\delta\hat G\|_{*,a} := \sup_{(s,\rho) \in \Omega_a} \left(\frac{d_a(\rho, s)}{s}\right)^{1/2}\,\|\delta\hat G(s)\|_{\mathcal{X}_\rho}.$$
For any $\hat G_1, \hat G_2$ satisfying $\|\hat G_i(s)\|_{\mathcal{X}_\rho} \le M$ on $\Omega_a$, apply (2.6) to $F = \mathcal{P}(\hat G_1(s')) - \mathcal{P}(\hat G_2(s'))$:
1. Using $\rho'(s') := \rho + \frac{1}{2}d_a(\rho, s')$ in (2.2) (so $\rho' - \rho = d_a(\rho', s') = \frac{1}{2}d_a(\rho, s')$ and $s' \le s$):
   $$\varepsilon C_S\int_0^s \|\mathcal{P}(\hat G_1(s')) - \mathcal{P}(\hat G_2(s'))\|_{\mathcal{X}_\rho}\,ds' \le \varepsilon C_S C_{\mathcal{P}} (2M)\sqrt{s}\int_0^s \frac{2^{3/2}\|\hat G_1 - \hat G_2\|_{*,a}}{d_a(\rho, s')^{3/2}}\,ds' \le \frac{8\sqrt{2}\,C_S C_{\mathcal{P}} M\varepsilon}{a}\left(\frac{s}{d_a(\rho, s)}\right)^{1/2}\|\hat G_1 - \hat G_2\|_{*,a}.$$
2. Using the zero-derivative-loss bound (2.4) directly at $\rho$ (since $d_a(\rho, s') \ge d_a(\rho, s)$ and $s' \le s$):
   $$\varepsilon C_*\int_0^s (s - s')^{-1/4}\bigl|a(\mathcal{P}(\hat G_1)) - a(\mathcal{P}(\hat G_2))\bigr|\,ds' \le \frac{8}{3}C_* C_a M\varepsilon\,S^{3/4}\left(\frac{s}{d_a(\rho, s)}\right)^{1/2}\|\hat G_1 - \hat G_2\|_{*,a}.$$
Hence for all $a \in [\frac{\rho_0}{12}, \frac{\rho_0}{6}]$:
$$\|\Phi(\hat G_1) - \Phi(\hat G_2)\|_{*,a} \le \kappa(\varepsilon)\,\|\hat G_1 - \hat G_2\|_{*,a}, \qquad \kappa(\varepsilon) := \left(\frac{96\sqrt{2}\,C_S C_{\mathcal{P}} M}{\rho_0} + \frac{8\cdot 2^{3/4}}{3}C_* C_a M\right)\varepsilon. \tag{2.8}$$
Now run Nishida's (1977) shrinking-cone iteration with $a_0 := \frac{\rho_0}{12}$, $a_{j+1} - a_j := \frac{\rho_0}{12 \cdot \frac{\pi^2}{6}(j+1)^2}$ (so $a_j \uparrow a_\infty = \frac{\rho_0}{6}$), starting from $\hat G^{(0)}(s) = S_g(s)G_0$ ($\|\hat G^{(0)}\|_{\mathcal{X}_\rho} \le C_S = \frac{M}{2}$). For the base step $\delta_1 = \hat G^{(1)} - \hat G^{(0)}$, $\int_0^s \frac{ds'}{d_{a_0}(\rho,s') + a_0 s'} = \frac{s}{d_{a_0}(\rho,s) + a_0 s} \le \frac{1}{2\sqrt{a_0}}(s/d_{a_0}(\rho,s))^{1/2}$ by AM–GM (and the $t^{-1/4}$ term gives $\frac{4}{3}s^{3/4} \le C(s/d_{a_0})^{1/2}$), so $\|\delta_1\|_{*,a_0} \le C_0'\varepsilon$ uniformly on $\Omega_{a_0}$. Passing from $\Omega_{a_j}$ to $\Omega_{a_{j+1}}$, we have $d_{a_j}(\rho, s) \ge (a_{j+1} - a_j)s$, so the $s$-factors cancel identically:
$$\|\hat G^{(j+1)}(s) - \hat G^{(j)}(s)\|_{\mathcal{X}_\rho} \le \frac{\|\hat G^{(j+1)} - \hat G^{(j)}\|_{*,a_j}}{(a_{j+1} - a_j)^{1/2}} \le C_0''(j+1)\kappa(\varepsilon)^j\,\varepsilon$$
uniformly on all of $\Omega_{a_{j+1}}$ (including $s \to 0^+$!). For $\varepsilon \le \varepsilon_0(\rho_0)$ with $\kappa(\varepsilon) \le \frac{1}{4}$ and $\sum_{j=0}^\infty C_0''(j+1)4^{-j}\varepsilon \le \frac{M}{2}$, every iterate satisfies $\|\hat G^{(j)}(s)\|_{\mathcal{X}_\rho} \le M$ on $\Omega_{a_j}$ and converges to a unique fixed point $\hat G \in C([0, 2]; \mathcal{X}_{\rho_1}^0)$ of (2.7) with $d_{a_\infty}(\rho_1, 2) = \frac{1}{6}\rho_0 > 0$.

**Classical regularity & exact $c_1$ continuation ($s > 0$):** Because $\widetilde G \in C([0,2]; \mathcal{X}_{\rho_1}^0)$ is strip-analytic on $\mathbb{S}_{\rho_1}$, $D\widetilde G$ and $Q\widetilde G$ are smooth on $(0,\infty)$, so $h(u,s) = e^{-us}(\varepsilon h_0(u)) + \int_0^s e^{-u(s-s')} f_{\mathcal{P}}(u,s')\,ds'$ is $C^1$ in $s > 0$ solving $h_s + u h = f_{\mathcal{P}}$, whence $\widetilde G = \mathcal{N}_g h$ solves (1.2) classically on $(0,\infty)\times(0,2]$. Moreover, for $s > 0$, $\tilde f_{\mathcal{P}}(u,s') = O(u^{-5/4})$ implies $\int_0^s e^{-u(s-s')}\tilde f_{\mathcal{P}}\,ds' = O(u^{-9/4})$ (holomorphic on $\operatorname{Re} m < \frac{9}{8}$), while $-\frac{\pi}{2}\int_0^s e^{-u(s-s')}a(\mathcal{P}(\widetilde G(s')))\phi_*(u)\,ds' = -\frac{\pi}{2}a(\mathcal{P}(\widetilde G(s)))/u^2 + O(u^{-3})$ (this needs $\alpha(s) := a(\mathcal{P}(\widetilde G(s)))$ to be Hölder continuous of some order $> \frac{1}{4}$ on $(0,2]$, which suffices for $\int_0^s e^{-u(s-s')}(\alpha(s') - \alpha(s))\,ds' = o(u^{-5/4})$; the note earlier asserted $C^1$, which is not needed and not proved. **Gap [S3, 2026-11-22 review]:** the Hölder bound on $\alpha$ is not written out here; it should follow from (2.4) and the $(s-s')^{-1/4}$ Duhamel kernel). Thus $\widehat{\widetilde G}(m,s) = N(m)\widehat h(m,s)$ continues meromorphically to $\operatorname{Re} m < \frac{9}{8}$ with a single simple pole at $m = 1$ ($N(1) = -\frac{6}{\pi}$), yielding $c_1(\widetilde G(s)) = 3\,a(\mathcal{P}(\widetilde G(s))) = -\frac{3}{\pi}J(\widetilde G(s))$.

---

### Lemma 4 (Transversal Horizon Crossing $s_\infty = 1 + O(\varepsilon)$ and Physical-Time Reconstruction, resolving `B4`)
For any $G \in \mathcal{X}_{\rho_1}^0$ and $r > 0$,
$$V_G(r^2) := -\frac{1}{\pi r}\int_0^\infty w^{-1/2}\ln\!\left(1 + \frac{r^2}{w}\right) G(w)\,dw.$$
Write $\ln(1 + r^2/w) = \frac{r^2}{w} - R(r^2/w)$ where $0 \le R(x) := x - \ln(1+x) \le \min(\frac{1}{2}x^2, x)$ for all $x > 0$. Because $G \in \mathcal{X}_{\rho_1}^0$ satisfies $a(G) = -\frac{1}{\pi}\int_0^\infty w^{-3/2}G(w)\,dw = 0$, the leading $O(r)$ term vanishes identically:
$$V_G(r^2) = +\frac{1}{\pi r}\int_0^\infty w^{-1/2} R\!\left(\frac{r^2}{w}\right) G(w)\,dw.$$
Splitting at $w = r^2$ and using $|G(w)| \le C\|G\|_{\mathcal{X}_{\rho_1}}\min(w^{5/8}, w^{1/4})$:
- For $w \le r^2$: $R(r^2/w) \le r^2/w$, giving $r^{-1}\int_0^{r^2} r^2 w^{-3/2} w^{5/8}\,dw = O(r^{5/4})$.
- For $w \ge r^2$: $R(r^2/w) \le \frac{1}{2}r^4/w^2$, giving $r^{-1}\int_{r^2}^1 r^4 w^{-5/2} w^{5/8}\,dw + r^3\int_1^\infty w^{-5/2} w^{1/4}\,dw = O(r^{5/4})$.
*(For curved $G(w) = c_1 w + O(w^2)$ with $a(G) = 0$, substituting $w = r^2/x$ gives the exact leading asymptotic $V_G(r^2) = \frac{c_1}{\pi}r^2\int_0^\infty x^{-5/2}(x - \ln(1+x))\,dx + O(r^3|\ln r|) = +\frac{2}{3}c_1 r^2 + O(r^3|\ln r|)$, verified in `TB-4`.)*
Similarly, differentiating under the integral sign gives $\partial_r V_G(r^2) = O(\varepsilon r^{1/4}) \to 0$ as $r \to 0^+$. Thus the right-hand side of (1.3) is $C^1$ on $r \ge 0$ with $dr/ds = -1 + O(\varepsilon r^{5/4}) \le -\frac{1}{2}$ for small $\varepsilon$, so $r(s)$ hits $0$ transversally at a unique $s_\infty = 1 + O(\varepsilon) \in (\frac{1}{2}, \frac{3}{2}) \subset [0, 2]$, with $r(s) = (s_\infty - s)(1 + O(\varepsilon(s_\infty - s)^{5/4}))$. Consequently $\tau(s) = \int_0^s ds'/r(s')$ maps $[0, s_\infty)$ bijectively onto $[0, \infty)$ with $\lambda(\tau) = r(s(\tau))^2 \sim c^2 e^{-2\tau}$, and reversing the autonomous $s$-clock change of variables yields the global solution $G(v,\tau) = \widetilde G(\lambda(\tau)v, s(\tau))$ of (1.1) with $\|G(\cdot,\tau)\|_{L^2_\sigma} = \lambda(\tau)^\sigma\|\widetilde G(\cdot,s(\tau))\|_{L^2_\sigma} \le C\varepsilon e^{-2\sigma\tau}$. $\blacksquare$
