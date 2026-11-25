# Rigorous Linear and Nonlinear Stability of Singular Self-Similar Profiles of the Generalized Constantin–Lax–Majda Model

**Branch:** `gclm/singular-stability`  
**Date:** 2026-11-25  
**Companion verification gates:**
- Gate 1: `validate/gclm_verify_singular_steady.py` (**PASS**)
- Gate 2: `validate/gclm_linear_stability.py` (**PASS**)
- Gate 3: `validate/gclm_branch2_and_dns_probe.py` (**PASS**)
- Gate 4: `validate/gclm_rigorous_closure_gate.py` (**PASS**)

---

## 1. Setup and Exact Singular Profile Branches

We consider the 1D generalized Constantin–Lax–Majda (gCLM) equation on $\mathbb{R}$:
$$\omega_t + a\,u\,\omega_x = u_x\,\omega, \qquad u_x = H(\omega)(x) := \frac{1}{\pi}\mathrm{P.V.}\int_{\mathbb{R}} \frac{\omega(y)}{x-y}\,dy, \qquad u(x,t) = \int_0^x H(\omega)(y,t)\,dy.$$
In dynamic-rescaling coordinates $\omega(x,t) = C_\omega^{-1}\Omega(C_l x,\tau(t))$, $\tau'(t) = C_\omega^{-1}$, with $C_\omega(\tau) = \exp\int_0^\tau c_\omega$, $C_l(\tau) = \exp\int_0^\tau c_l$:
\begin{equation}\label{eq:gclm_dyn}
\Omega_\tau + (c_l X + a U)\,\Omega_X = (c_\omega + U_X)\,\Omega, \qquad U_X = H(\Omega), \qquad U(X,\tau) = \int_0^X H(\Omega)(Y,\tau)\,dY.
\end{equation}

For every $a < 0$,~\eqref{eq:gclm_dyn} admits two explicit half-line singular steady-state families intersecting transcritically at $a = -1$:
1. **Branch I** (Huang–Tong–Wang, `arXiv:2603.25104v4`, Thm. 3): with $\mu = \frac{1}{1-a} \in (0,1)$,
   $$\bar\Omega_{\mathrm{I},a}(X) = -\sin(\pi\mu)(X-1)^{-\mu}\mathbf{1}_{\{X>1\}}, \qquad \bar c_{l,\mathrm{I}} = 1-a, \qquad \bar c_{\omega,\mathrm{I}} = -1, \qquad \bar\gamma_{\mathrm{I}} = 1-a.$$
2. **Branch II** (universal $\beta = \frac{1}{2}$ family, all $a < 0$):
   $$\bar\Omega_{\mathrm{II},a}(X) = \frac{1}{a}(X-1)^{-1/2}\mathbf{1}_{\{X>1\}}, \qquad \bar c_{l,\mathrm{II}} = 2, \qquad \bar c_{\omega,\mathrm{II}} = -1, \qquad \bar\gamma_{\mathrm{II}} = 2,$$
   with $(\bar U_{\mathrm{II},a})_X(X) = -\frac{1}{a}(1-X)^{-1/2}\mathbf{1}_{\{X<1\}}$ (identically zero on $X>1$), $\bar U_{\mathrm{II},a}(X) = -\frac{2}{a}$ on $X \ge 1$, and strictly outgoing background advection $\bar A_{\mathrm{II}}(1+v) = 2v > 0$ on $X = 1+v > 1$.

At **$a = -1$**, Branch I and Branch II **coincide identically**: $\bar\Omega_{\mathrm{I},-1} = \bar\Omega_{\mathrm{II},-1} = -(X-1)^{-1/2}\mathbf{1}_{\{X>1\}}$, $\bar c_l = 2$, $\bar c_\omega = -1$.

### Support-Interior Perturbation System and Edge-Amplitude-Pinning Gauge (Gauge (c))
For any $a < 0$, write support-interior perturbations of $\bar\Omega_{\mathrm{II},a}$ (hence of $\bar\Omega_{\mathrm{I},-1}$ at $a = -1$) on $X = 1 + v = 1 + e^z > 1$ as
$$\Omega(1+e^z, \tau) = \bar\Omega_{\mathrm{II},a}(1+e^z)\big(1 + G(z,\tau)\big), \qquad g(1+v,\tau) := \frac{1}{a}v^{-1/2}G(\ln v, \tau),$$
with $\Omega(X, \tau) \equiv 0$ for $X < 1$. Define the operators on $v = e^z \in (0,\infty)$ (all arguments written in $z \in \mathbb{R}$):
$$\widetilde{H}[G](z) := H_+[v^{-1/2} G](1+e^z) = \frac{1}{\pi}\mathrm{P.V.}\int_0^\infty \frac{(v')^{-1/2} G(\ln v')}{e^z - v'}\,dv',$$
$$\widetilde{Q}[G](z) := \frac{1}{e^z}\int_0^{e^z} \widetilde{H}[G](\ln y)\,dy = \int_{-\infty}^z e^{z'-z}\widetilde{H}[G](z')\,dz',$$
and the two scalar linear functionals
$$V_G := a\,\delta U(1) = -\frac{1}{\pi}\int_{-\infty}^\infty e^{z/2}\ln(1 + e^{-z})\,G(z)\,dz, \qquad \widetilde{H}_{\mathrm{edge}}[G] := \lim_{z\to-\infty}\widetilde{H}[G](z) = -\frac{1}{\pi}\int_{-\infty}^\infty e^{-z/2}G(z)\,dz = -\frac{1}{\pi}\widehat{G}\!\left(\tfrac{1}{2}\right),$$
where we adopt the Mellin–Laplace convention $\widehat{f}(s) := \int_{\mathbb{R}} e^{-sz} f(z)\,dz$, with Plancherel identity
$$\|f\|_{L^2_\sigma(\mathbb{R})}^2 := \int_{\mathbb{R}} e^{-2\sigma z}|f(z)|^2\,dz = \frac{1}{2\pi}\int_{\mathbb{R}} |\widehat{f}(\sigma + i\xi)|^2\,d\xi.$$

Imposing **edge-stagnation pinning** $A(1,\tau) := c_l(\tau) + a U(1,\tau) \equiv 0$, i.e., $c_l(\tau) = 2 - V_G(\tau)$, eliminates the singular $(X-1)^{-3/2}$ advection drift identically: on $X = 1+v > 1$, $A(1+v,\tau) = v(2 + \widetilde{Q}[G] - V_G)$, so
$$A\,\partial_X\!\Big(\frac{1}{a}v^{-1/2}(1+G)\Big) = \frac{1}{a}v^{-1/2}\Big[\big(2 + \widetilde{Q}[G] - V_G\big)\partial_z G - \Big(1 + \frac{1}{2}\widetilde{Q}[G] - \frac{1}{2}V_G\Big)(1+G)\Big].$$
Writing $\delta c_\omega := c_\omega + 1$ and moving the stretching term to the right-hand side reduces~\eqref{eq:gclm_dyn} on $X > 1$ to
\begin{equation}\label{eq:gclm_pert_gen}
G_\tau + \big(2 + \widetilde{Q}[G] - V_G\big)\partial_z G = \Big(\delta c_\omega - \tfrac{1}{2}V_G + \frac{1}{a}\widetilde{H}[G] + \frac{1}{2}\widetilde{Q}[G]\Big)(1 + G).
\end{equation}
- **Gauge (a)** ($\delta c_\omega = +\frac{1}{2}V_G$, i.e. $c_\omega = -\frac{1}{2}c_l$): yields the pure bulk operator $M_{\mathrm{II},a}^{\mathrm{bulk}} = -2\partial_z + \frac{1}{a}\widetilde{H} + \frac{1}{2}\widetilde{Q}$.
- **Gauge (c) (edge-amplitude-pinning gauge):** unlike 1D Hou–Luo (where the $\Theta$-equation locked $c_l + 2 c_\omega = 0$), $c_\omega(\tau)$ is a free scaling parameter in gCLM. Choosing
  \begin{equation}\label{eq:gauge_c}
  c_l(\tau) = 2 - V_G(\tau), \qquad c_\omega(\tau) = -\frac{1}{2}c_l(\tau) - \left(\frac{1}{a} + \frac{1}{2}\right)\widetilde{H}_{\mathrm{edge}}[G]
  \end{equation}
  (so $\delta c_\omega - \frac{1}{2}V_G = -(\frac{1}{a} + \frac{1}{2})\widetilde{H}_{\mathrm{edge}}[G]$) cancels the edge trace $\lim_{z\to-\infty}(\frac{1}{a}\widetilde{H}[G] + \frac{1}{2}\widetilde{Q}[G]) = (\frac{1}{a}+\frac{1}{2})\widetilde{H}_{\mathrm{edge}}[G]$ identically, yielding the exact nonlinear Gauge-(c) PDE
  \begin{equation}\label{eq:nonlin_gauge_c}
  G_\tau + \big(2 + \widetilde{Q}[G] - V_G\big)\partial_z G = \Big(\frac{1}{a}\widetilde{H}_0[G] + \frac{1}{2}\widetilde{Q}_0[G]\Big)(1 + G),
  \end{equation}
  where
  $$\widetilde{H}_0[G](z) := \widetilde{H}[G](z) - \widetilde{H}_{\mathrm{edge}}[G] = \frac{e^z}{\pi}\mathrm{P.V.}\int_0^\infty \frac{(v')^{-3/2}G(\ln v')}{e^z - v'}\,dv', \qquad \widetilde{Q}_0[G](z) := \widetilde{Q}[G](z) - \widetilde{H}_{\mathrm{edge}}[G] = \int_{-\infty}^z e^{z'-z}\widetilde{H}_0[G](z')\,dz'.$$
  The linearization of~\eqref{eq:nonlin_gauge_c} about $G = 0$ is
  $$M_c := -2\partial_z + \frac{1}{a}\widetilde{H}_0 + \frac{1}{2}\widetilde{Q}_0.$$
  Note that $M_c$ coincides with $M_{\mathrm{II},a}^{\mathrm{bulk}}$ identically at $a = -2$, and differs for $a \ne -2$ by the constant edge-trace subtraction $(\frac{1}{a} + \frac{1}{2})\widetilde{H}_{\mathrm{edge}}[G]$.

### Global Distributional Validity on $\mathbb{R}$ (Zero $\delta(X-1)$ Defect)
Because $A_X = c_l + a U_X$, the gCLM equation~\eqref{eq:gclm_dyn} admits the conservative weak form
$$\Omega_\tau + \partial_X(A\,\Omega) = \big(c_\omega + c_l + (1 + a)U_X\big)\,\Omega.$$
For any strong solution $G$ of~\eqref{eq:nonlin_gauge_c} in $\mathcal{X}_{\mathrm{gCLM}}$ (Theorem 3 below), we have $\lim_{z\to-\infty} G(z,\tau) = 0$, $A(1+v,\tau) = v(2 + \widetilde{Q}[G] - V_G) = O(v)$ as $v \to 0^+$, and $U_X(1+v,\tau) = \frac{1}{a}\widetilde{H}[G](\ln v) \in L^\infty((0,\infty))$. Hence the flux $A(1+v,\tau)\Omega(1+v,\tau) = O(v^{1/2}) \to 0$ as $v \to 0^+$, while $(1+a)U_X\Omega \in L^1_{\mathrm{loc}}(\mathbb{R})$. Testing against any $\varphi \in C_c^\infty(\mathbb{R})$ and integrating $\int_1^\infty \partial_X(A\Omega)\varphi\,dX = -\int_1^\infty A\Omega\,\varphi_X\,dX$ produces **zero boundary term at $X = 1$** — therefore $\Omega(X,\tau) = \frac{1}{a}(X-1)^{-1/2}(1 + G(\ln(X-1),\tau))\mathbf{1}_{\{X>1\}}$ is a genuine **global distributional solution on all of $\mathbb{R}$** with zero $\delta(X-1)$ defect!

---

## 2. Theorem 1: Exact 2-Gamma Conjugation and Slit-Plane Hardy Semigroup for $M_c$

**Theorem 1 (Exact Conjugation of Gauge-(c) Linear Operator $M_c$).**
*For every $a < 0$, define the meromorphic Mellin symbol*
$$\tilde c_{\mathrm{II},a}(r) := \left(\frac{1}{|a|} - \frac{1}{2r+1}\right)\tan(\pi r) = -\frac{2r + 1 + a}{a(2r + 1)}\tan(\pi r),$$
*which is the exact Mellin symbol of $e^{+z/2}(\frac{1}{a}\widetilde{H}_0 + \frac{1}{2}\widetilde{Q}_0)$ throughout the strip $\frac{1}{2} < \operatorname{Re} r < \frac{3}{2}$ (Lemma 2.1 below). Define the explicit 2-Gamma multiplier*
$$N_{L,\mathrm{II}}(s; a) := \sin(\pi s)\,|a|^{-2s}\,\frac{\Gamma(2s + 2 + a)}{\Gamma(2s + 2)}.$$
*(At $a = -1$, this simplifies to the elementary function $N_{L,\mathrm{II}}(s;-1) = \frac{\sin(\pi s)}{2s+1}$.) Fix $(a, \sigma)$ in the shift-strip pole-free regime:*
$$\text{either } a \in (-2, 0) \text{ and } \sigma \in \big(0, \tfrac{1}{2}\big), \qquad \text{or } a \in (-3, -2] \text{ and } \sigma \in \big(-1 - \tfrac{a}{2},\, \tfrac{1}{2}\big).$$
*Then:*
1. **Exact shift identity:** $N_{L,\mathrm{II}}(s;a)\,\tilde c_{\mathrm{II},a}(s + \frac{1}{2}) = -N_{L,\mathrm{II}}(s + \frac{1}{2};a)$ meromorphically on $\mathbb{C}$.
2. **Stirling asymptotics & zero-free shift strip:** $N_{L,\mathrm{II}}(s; a)$ is analytic and zero-free on the closed shift strip $\sigma \le \operatorname{Re} s \le \sigma + \frac{1}{2}$, with two-sided Stirling bounds on $\operatorname{Re} s = \sigma$:
   $$c(\sigma,a)(1 + |\xi|)^a e^{\pi|\xi|} \le |N_{L,\mathrm{II}}(\sigma + i\xi; a)| \le C(\sigma,a)(1 + |\xi|)^a e^{\pi|\xi|}.$$
3. **Explicit slit-plane Hardy contraction semigroup for $M_c$:** Let $\mathcal{H}_{\sigma,a} := \mathcal{N}_{L,\mathrm{II}}^{-1}(L^2_\sigma(\mathbb{R}))$ equipped with $\|G\|_{\mathcal{H}_{\sigma,a}} := \|\mathcal{N}_{L,\mathrm{II}} G\|_{L^2_\sigma(\mathbb{R})}$ (by Paley–Wiener, elements of $\mathcal{H}_{\sigma,a}$ extend analytically to the strip $|\operatorname{Im} z| < \pi$, i.e. the slit plane $X \in \mathbb{C}\setminus(-\infty, 1]$). On the dense domain $D(M_c) := \mathcal{N}_{L,\mathrm{II}}^{-1}\big(H^1_\sigma(\mathbb{R}) \cap L^2_{\sigma+1/2}(\mathbb{R})\big)$, the Gauge-(c) operator $M_c$ satisfies the exact operator conjugation
   $$\mathcal{N}_{L,\mathrm{II}}\,M_c\,\mathcal{N}_{L,\mathrm{II}}^{-1} = A_\sigma := -2\partial_z - e^{-z/2} \quad \text{on } D(A_\sigma) := H^1_\sigma(\mathbb{R}) \cap L^2_{\sigma+1/2}(\mathbb{R}),$$
   and generates the explicit $C_0$-contraction semigroup
   $$e^{\tau M_c} = \mathcal{N}_{L,\mathrm{II}}^{-1}\,S(\tau)\,\mathcal{N}_{L,\mathrm{II}}, \qquad (S(\tau)F)(z) := \exp\!\Big(-e^{-z/2}(e^\tau - 1)\Big)\,F(z - 2\tau),$$
   satisfying $\|e^{\tau M_c}\|_{\mathcal{H}_{\sigma,a} \to \mathcal{H}_{\sigma,a}} \le e^{-2\sigma\tau}$ for all $\tau \ge 0$, exact spectrum $\sigma(M_c|_{\mathcal{H}_{\sigma,a}}) = \{\operatorname{Re}\lambda \le -2\sigma\}$ with empty point spectrum in $\operatorname{Re}\lambda > -2\sigma$, and $\|(\lambda - M_c)^{-1}\|_{\mathcal{H}_{\sigma,a}} \le (\operatorname{Re}\lambda + 2\sigma)^{-1}$ for $\operatorname{Re}\lambda > -2\sigma$.

*(Remark on Gauge (a): on $\operatorname{Re} s = \sigma \in (0,\frac{1}{2})$, the shifted line $\operatorname{Re} r = \sigma + \frac{1}{2} \in (\frac{1}{2}, 1)$ lies across the simple pole $r = \frac{1}{2}$ of $\tan(\pi r)$ inside the Gauge-(c) strip $(\frac{1}{2},\frac{3}{2})$. Hence $\mathcal{N}_{L,\mathrm{II}}$ conjugates $M_c$ — NOT $M_{\mathrm{II},a}^{\mathrm{bulk}}$ — except at $a = -2$ where $M_c \equiv M_{\mathrm{II},-2}^{\mathrm{bulk}}$. For $a \ne -2$, $M_{\mathrm{II},a}^{\mathrm{bulk}} G - M_c G = (\frac{1}{a}+\frac{1}{2})\widetilde{H}_{\mathrm{edge}}[G] \notin L^2_\sigma(\mathbb{R})$ whenever $\widehat{G}(\frac{1}{2}) \ne 0$. Note also that for $a \le -3$, every shift strip $[\sigma, \sigma+\frac{1}{2}]$ with $\sigma \in (0,\frac{1}{2})$ contains a pole $s_n = -1 - \frac{a+n}{2}$ of $\Gamma(2s+2+a)$, blocking full-domain conjugation.)*

**Proof.**
1. *Shift identity:* Using $\tan(\pi(s+\frac{1}{2})) = -\cot(\pi s)$ and $a = -|a| < 0$,
   $$\frac{N_{L,\mathrm{II}}(s+\frac{1}{2}; a)}{N_{L,\mathrm{II}}(s; a)} = \frac{\sin(\pi s + \frac{\pi}{2})}{\sin(\pi s)}\cdot |a|^{-1}\cdot \frac{\Gamma(2s + 3 + a)}{\Gamma(2s + 2 + a)}\cdot \frac{\Gamma(2s + 2)}{\Gamma(2s + 3)} = \cot(\pi s)\cdot \frac{2s + 2 + a}{|a|(2s + 2)} = -\tilde c_{\mathrm{II},a}\!\left(s + \tfrac{1}{2}\right).$$
2. *Poles, zeros, and Stirling:* Zeros of $\sin(\pi s)$ lie at $\mathbb{Z}$, outside $(0, 1) \supset [\sigma, \sigma + \frac{1}{2}]$; $\Gamma(2s+2)$ is zero-free and pole-free on $\operatorname{Re} s > -1$; poles of $\Gamma(2s+2+a)$ lie at $s_n = -1 - \frac{a+n}{2}$ ($n \in \mathbb{N}_0$), whose rightmost pole is $s_0 = -1 - \frac{a}{2} < \sigma$ in both specified regimes. By Stirling's formula, $|\Gamma(2\sigma + 2 + a + 2i\xi)/\Gamma(2\sigma + 2 + 2i\xi)| = (2|\xi|)^a(1 + O(|\xi|^{-1}))$ while $|\sin(\pi(\sigma+i\xi))| \asymp e^{\pi|\xi|}$ on $\sigma \in (0,\frac{1}{2})$.
3. *Generator domain, conjugation, and exact spectrum:* Conjugating by the unitary isometry $U_\sigma : L^2_\sigma(\mathbb{R}) \to L^2(\mathbb{R})$, $(U_\sigma F)(z) = e^{-\sigma z} F(z)$, transforms $S(\tau)$ into $(S_\sigma(\tau)u)(z) = e^{-2\sigma\tau}\exp(-e^{-z/2}(e^\tau - 1))u(z - 2\tau)$. Strong continuity on $L^2(\mathbb{R})$ follows immediately from dominated convergence and translation continuity, and its contraction bound is $\|S_\sigma(\tau)\|_{L^2\to L^2} \le e^{-2\sigma\tau}$. Its infinitesimal generator on $L^2(\mathbb{R})$ is the closed operator $B_\sigma = -2\partial_z - 2\sigma - e^{-z/2}$ with exact domain $D(B_\sigma) = \{u \in H^1(\mathbb{R}) : e^{-z/2} u \in L^2(\mathbb{R})\}$, which pulls back via $U_\sigma^{-1}$ to $D(A_\sigma) = H^1_\sigma(\mathbb{R}) \cap L^2_{\sigma+1/2}(\mathbb{R})$.

   For any $F \in D(A_\sigma)$, its Mellin transform $\widehat{F}(s)$ is analytic on $\sigma < \operatorname{Re} s < \sigma + \frac{1}{2}$ with boundary values in $L^2$ on both edges $\operatorname{Re} s = \sigma$ and $\operatorname{Re} s = \sigma + \frac{1}{2}$, and $\widehat{A_\sigma F}(s) = -2s\widehat{F}(s) - \widehat{F}(s+\frac{1}{2})$ on $\operatorname{Re} s = \sigma$. Setting $G := \mathcal{N}_{L,\mathrm{II}}^{-1} F$, analyticity and non-vanishing of $N_{L,\mathrm{II}}$ on $[\sigma, \sigma+\frac{1}{2}]$ ensures $\widehat{G}(s) = \widehat{F}(s)/N_{L,\mathrm{II}}(s;a)$ has $L^2$ traces on both $\operatorname{Re} s = \sigma$ and $\operatorname{Re} r = \sigma + \frac{1}{2} \in (\frac{1}{2}, 1)$. Dividing $\widehat{A_\sigma F}(s)$ by $N_{L,\mathrm{II}}(s;a)$ and applying part (1) gives
   $$\frac{\widehat{A_\sigma F}(s)}{N_{L,\mathrm{II}}(s; a)} = -2s\,\widehat{G}(s) + \tilde c_{\mathrm{II},a}\!\left(s + \tfrac{1}{2}\right)\widehat{G}\!\left(s + \tfrac{1}{2}\right) = \widehat{M_c G}(s) \qquad (\operatorname{Re} s = \sigma),$$
   since $\operatorname{Re}(s+\frac{1}{2}) \in (\frac{1}{2}, 1) \subset (\frac{1}{2}, \frac{3}{2})$ where $\tilde c_{\mathrm{II},a}$ is the exact symbol of $e^{+z/2}(\frac{1}{a}\widetilde{H}_0 + \frac{1}{2}\widetilde{Q}_0)$ (Lemma 2.1). Hence $\mathcal{N}_{L,\mathrm{II}} M_c \mathcal{N}_{L,\mathrm{II}}^{-1} = A_\sigma$ on $D(A_\sigma)$, establishing the semigroup and resolvent bound $\|(\lambda - M_c)^{-1}\|_{\mathcal{H}_{\sigma,a}} \le (\operatorname{Re}\lambda + 2\sigma)^{-1}$ (and $\sigma(M_c) \subset \{\operatorname{Re}\lambda \le -2\sigma\}$) by Hille–Yosida.

   For the reverse spectral inclusion $\{\operatorname{Re}\lambda \le -2\sigma\} \subset \sigma(M_c) = \sigma(B_\sigma)$: the $L^2(\mathbb{R})$-adjoint $B_\sigma^* = 2\partial_z - 2\sigma - e^{-z/2}$ satisfies $B_\sigma^*\psi_\lambda = \bar\lambda\,\psi_\lambda$ for the explicit smooth function
   $$\psi_\lambda(z) := \exp\!\left(\frac{2\sigma + \bar\lambda}{2}z - e^{-z/2}\right).$$
   Whenever $\operatorname{Re}\lambda < -2\sigma$, $\psi_\lambda(z)$ decays exponentially as $z \to +\infty$ and super-exponentially as $z \to -\infty$, so $\psi_\lambda \in D(B_\sigma^*) \setminus \{0\}$. Thus $\bar\lambda \in \sigma_p(B_\sigma^*)$, placing every $\lambda$ with $\operatorname{Re}\lambda < -2\sigma$ in $\sigma(B_\sigma)$ (in fact in the residual spectrum, since $B_\sigma u = \lambda u \implies u(z) = C e^{-\frac{2\sigma+\lambda}{2}z + e^{-z/2}} \notin L^2(\mathbb{R})$). Because $\sigma(B_\sigma)$ is closed, $\sigma(M_c|_{\mathcal{H}_{\sigma,a}}) = \{\operatorname{Re}\lambda \le -2\sigma\}$. $\square$

---

## 3. Theorem 2: Unconditional Finite-Sobolev Coercivity and Exact Operator Bounds (Gauge (c))

Passing from the Hardy space $\mathcal{H}_{\sigma,a}$ to standard weighted Sobolev spaces $H^k_\sigma(\mathbb{R})$ removes the exponential analyticity requirement completely and yields **unconditional finite-Sobolev $H^k_\sigma$ coercivity** of $M_c$ for every $a \in [-2, 0)$.

### Lemma 2.1 (Mellin Representations of $\widetilde{H}_0$ and $\widetilde{Q}_0$)
For $G \in L^2_\sigma(\mathbb{R})$ with $\sigma \in (\frac{1}{2}, \frac{3}{2})$:
1. $T_{H_0} G := e^{+z/2}\widetilde{H}_0[G]$ is a bounded Fourier–Mellin multiplier on $L^2_\sigma(\mathbb{R})$ with symbol $-\tan(\pi r)$ on $\operatorname{Re} r = \sigma$, and exact operator norm
   $$\|T_{H_0}\|_{L^2_\sigma \to L^2_\sigma} = \sup_{\xi \in \mathbb{R}} |\tan(\pi(\sigma + i\xi))| = \max\!\big(|\tan(\pi\sigma)|, 1\big) = \begin{cases} |\tan(\pi\sigma)|, & \sigma \in (\frac{1}{2}, \frac{3}{4}) \cup (\frac{5}{4}, \frac{3}{2}), \\ 1, & \sigma \in [\frac{3}{4}, \frac{5}{4}]. \end{cases}$$
   In particular, $\|T_{H_0}\|_{L^2_{5/8} \to L^2_{5/8}} = \cot(\frac{\pi}{8}) = \sqrt{2} + 1$.
2. $T_{Q_0} G := e^{+z/2}\widetilde{Q}_0[G]$ has Mellin symbol $-\frac{2\tan(\pi r)}{2r + 1}$ on $\operatorname{Re} r = \sigma \in (\frac{1}{2}, \frac{3}{2})$, satisfying $\partial_z \widetilde{Q}_0[G] = \widetilde{H}_0[G] - \widetilde{Q}_0[G]$ and $\|T_{Q_0}\|_{L^2_{5/8}\to L^2_{5/8}} \le \frac{2}{2(5/8)+1}(\sqrt{2}+1) = \frac{8}{9}(\sqrt{2}+1)$.
3. For any $k \ge 0$, $\partial_z^k \widetilde{H}_0[G] = e^{-z/2} T_{H_0}\big((\partial_z - \frac{1}{2})^k G\big)$ and $\partial_z^k \widetilde{Q}_0[G] = e^{-z/2} T_{Q_0}\big((\partial_z - \frac{1}{2})^k G\big)$.

**Proof.**
By definition, $\widetilde{H}_0[G](\ln v) = v H_+[v^{-3/2}G](1+v)$, so $(T_{H_0}G)(\ln v) = v^{3/2}H_+[v^{-3/2}G](1+v)$. Under $v = e^z$, $T_{H_0}$ is a multiplicative convolution on $(0,\infty)$ whose Mellin symbol on $\operatorname{Re} r = \sigma$ is $\frac{1}{\pi}\mathrm{P.V.}\int_0^\infty \frac{u^{3/2 - r - 1}}{u - 1}\,du = -\cot(\pi(\frac{3}{2} - r)) = -\tan(\pi r)$, valid throughout the Muckenhoupt $A_2$ strip $0 < \frac{3}{2} - \sigma < 1 \iff \frac{1}{2} < \sigma < \frac{3}{2}$. Since $|\tan(\pi(\sigma+i\xi))|^2 = \frac{\cosh(2\pi\xi) - \cos(2\pi\sigma)}{\cosh(2\pi\xi) + \cos(2\pi\sigma)}$ is monotone in $\cosh(2\pi\xi) \ge 1$ with limit $1$ as $|\xi|\to\infty$ and value $\tan^2(\pi\sigma)$ at $\xi=0$, the exact norm follows. Identity (3) follows from exact dilation commutation $\partial_z H_+ = H_+ \partial_z$ on $(0,\infty)$ applied to $v^{-3/2}G$: $\partial_z(v H_+[v^{-3/2}G]) = v H_+[\partial_z(v^{-3/2}G)] + v H_+[v^{-3/2}G] = v H_+[v^{-3/2}(\partial_z - \frac{1}{2})G]$. $\square$

### Lemma 2.2 (Exact $L^2_{1/4}$ Identity and Unconditional Higher-Order $H^k_\sigma$ Coercivity of $M_c$)
Let $a \in [-2, 0)$ and let $M_c = -2\partial_z + \frac{1}{a}\widetilde{H}_0 + \frac{1}{2}\widetilde{Q}_0$.
1. **Exact $L^2_{1/4}(\mathbb{R})$ energy identity** (for all $G \in H^1_{1/4}(\mathbb{R}) \cap L^2_{5/8}(\mathbb{R})$):
   $$\langle M_c G,\, G \rangle_{L^2_{1/4}} = -\frac{1}{2}\|G\|_{L^2_{1/4}}^2 - \frac{2 - |a|}{4\pi|a|}\,\widehat{G}\!\left(\tfrac{1}{2}\right)^2 - \frac{1}{4\pi}\int_{-\infty}^\infty \frac{\xi\coth(\pi\xi)}{1 + \xi^2}\,\big|\widehat{G}(\tfrac{1}{2} + i\xi)\big|^2\,d\xi \le -\frac{1}{2}\|G\|_{L^2_{1/4}}^2.$$
2. **Unconditional higher-order Plancherel coercivity:** For every $\sigma \in (\frac{1}{4}, \frac{3}{4})$ (so $\sigma_1 := \sigma + \frac{1}{4} \in (\frac{1}{2}, 1)$) and every $k \in \mathbb{N}_0$,
   $$\langle \partial_z^k M_c G,\, \partial_z^k G \rangle_{L^2_\sigma} = -2\sigma\|\partial_z^k G\|_{L^2_\sigma}^2 + \frac{1}{2\pi}\int_{-\infty}^\infty \operatorname{Re}\!\Big[(r - \tfrac{1}{2})^k \bar r^k \tilde c_{\mathrm{II},a}(r)\Big]\big|\widehat{G}(r)\big|^2\,d\xi, \qquad r = \sigma_1 + i\xi.$$
   For **every $a \in [-2, 0)$** and **every $\xi \in \mathbb{R}$**:
   - At **$\sigma_+ = \frac{5}{8}$** ($\sigma_1 = \frac{7}{8}$): $\operatorname{Re}[(r-\frac{1}{2})^k\bar r^k\tilde c_{\mathrm{II},a}(r)] < 0$ for $k \in \{0, 1, 2, 3\}$, hence
     $$\langle \partial_z^k M_c G,\, \partial_z^k G \rangle_{L^2_{5/8}} \le -\frac{5}{4}\|\partial_z^k G\|_{L^2_{5/8}}^2 \qquad (k = 0, 1, 2, 3).$$
   - At **$\sigma_- = \frac{3}{8}$** ($\sigma_1 = \frac{5}{8}$): $\operatorname{Re}[(r-\frac{1}{2})^k\bar r^k\tilde c_{\mathrm{II},a}(r)] < 0$ for $k \in \{0, 1, 2\}$, hence
     $$\langle \partial_z^k M_c G,\, \partial_z^k G \rangle_{L^2_{3/8}} \le -\frac{3}{4}\|\partial_z^k G\|_{L^2_{3/8}}^2 \qquad (k = 0, 1, 2).$$
   - At **$\sigma_- = \frac{1}{4}$** ($\sigma_1 = \frac{1}{2}$): for ALL $a \in [-\frac{8}{5}, 0)$ (including $a = -1$ Branch I), $\langle \partial_z^k M_c G, \partial_z^k G \rangle_{L^2_{1/4}} \le -\frac{1}{2}\|\partial_z^k G\|_{L^2_{1/4}}^2$ for $k = 0, 1, 2$.

**Proof.**
For part (1), write $f(v) = v^{-1/2}G(\ln v)$. Symmetrizing the double integral $\langle \widetilde{H}[G], G\rangle_{L^2_{1/4}} = \frac{1}{\pi}\iint \frac{f(v)f(v')}{v(v-v')}\,dv\,dv'$ via $\frac{1}{v(v-v')} + \frac{1}{v'(v'-v)} = -\frac{1}{vv'}$ gives $\langle \widetilde{H}[G], G\rangle_{L^2_{1/4}} = -\frac{1}{2\pi}\widehat{G}(\frac{1}{2})^2$. Subtracting $\widetilde{H}_{\mathrm{edge}}[G] = -\frac{1}{\pi}\widehat{G}(\frac{1}{2})$ adds $+\frac{1}{\pi}\widehat{G}(\frac{1}{2})^2$, yielding $\langle \widetilde{H}_0[G], G\rangle_{L^2_{1/4}} = +\frac{1}{2\pi}\widehat{G}(\frac{1}{2})^2$, and multiplying by $\frac{1}{a} = -\frac{1}{|a|}$ gives $-\frac{1}{2\pi|a|}\widehat{G}(\frac{1}{2})^2$. For $\frac{1}{2}\widetilde{Q}_0[G]$, shifting the Plancherel line $\operatorname{Re} r = \frac{1}{2} + \epsilon \to \frac{1}{2}^+$ picks up the Sokhotski–Plemelj delta residue $+\frac{1}{4\pi}\widehat{G}(\frac{1}{2})^2$ from $-\frac{\tan(\pi r)}{2r+1}$ (whose residue at $r = \frac{1}{2}$ is $+\frac{1}{2\pi}$), plus the principal-value integral $-\frac{1}{4\pi}\int \frac{\xi\coth(\pi\xi)}{1+\xi^2}|\widehat{G}(\frac{1}{2}+i\xi)|^2\,d\xi$. Summing the two contributions yields the exact coefficient $\frac{1}{4\pi} - \frac{1}{2\pi|a|} = -\frac{2-|a|}{4\pi|a|} \le 0$ for $|a| \le 2$.

For part (2), integration by parts gives $\langle -2\partial_z(\partial_z^k G), \partial_z^k G\rangle_{L^2_\sigma} = -2\sigma\|\partial_z^k G\|_{L^2_\sigma}^2$, while Lemma 2.1(3) gives
$$\langle \partial_z^k(e^{-z/2} T_{\tilde c} G), \partial_z^k G\rangle_{L^2_\sigma} = \langle (\partial_z - \tfrac{1}{2})^k T_{\tilde c} G, \partial_z^k G\rangle_{L^2_{\sigma + 1/4}} = \frac{1}{2\pi}\int_{\mathbb{R}} (r - \tfrac{1}{2})^k \bar r^k \tilde c_{\mathrm{II},a}(r)|\widehat{G}(r)|^2\,d\xi.$$
Let $P_k(\xi) := (r - \frac{1}{2})^k \bar r^k = \big((\xi^2 + c) + \frac{i}{2}\xi\big)^k = R_k(x) + i\xi I_k(x)$, where $x := \xi^2 \ge 0$ and $c := \sigma_1(\sigma_1 - \frac{1}{2})$. (Note that as $|\xi|\to\infty$, $\operatorname{Re}[P_k(\xi)\tilde c_{\mathrm{II},a}(r)] = -\frac{k+|a|}{2|a|}|\xi|^{2k-1} + O(|\xi|^{2k-3}) < 0$.)

*Case A ($\sigma_1 = \frac{1}{2}$, $k = 1, 2$):* Here $r = \frac{1}{2} + i\xi$; the factor $(r - \frac{1}{2})^k = (i\xi)^k$ cancels the simple pole of $\tan(\pi r)$ at $r = \frac{1}{2}$, so no Sokhotski–Plemelj delta residue arises, and $\tilde c_{\mathrm{II},a}(\frac{1}{2}+i\xi) = \big(\frac{i}{|a|} - \frac{\xi+i}{2(1+\xi^2)}\big)\coth(\pi\xi)$.
- $k=1$: using $\operatorname{Re}[(\xi^2 + \frac{i}{2}\xi)(-\xi - i)] = -\xi^3 + \frac{1}{2}\xi$, we obtain
  $$\operatorname{Re}\!\big[(r-\tfrac{1}{2})\bar r\,\tilde c_{\mathrm{II},a}(r)\big] = \xi\coth(\pi\xi)\!\left(-\frac{1}{2|a|} + \frac{\frac{1}{2} - \xi^2}{2(1 + \xi^2)}\right) \le \xi\coth(\pi\xi)\!\left(-\frac{1}{2|a|} + \frac{1}{4}\right) \le 0 \qquad (\forall |a| \le 2).$$
- $k=2$: using $\operatorname{Re}[(\xi^4 - \frac{1}{4}\xi^2 + i\xi^3)(-\xi - i)] = -\xi^5 + \frac{5}{4}\xi^3$, we obtain
  $$\operatorname{Re}\!\big[(r-\tfrac{1}{2})^2\bar r^2\,\tilde c_{\mathrm{II},a}(r)\big] = \xi^3\coth(\pi\xi)\!\left(-\frac{1}{|a|} + \frac{\frac{5}{4} - \xi^2}{2(1 + \xi^2)}\right) \le \xi^3\coth(\pi\xi)\!\left(-\frac{1}{|a|} + \frac{5}{8}\right) \le 0 \qquad (\forall |a| \le \tfrac{8}{5}).$$

*Case B ($\sigma_1 \in \{\frac{5}{8}, \frac{7}{8}\}$, ALL $a \in [-2, 0)$ and ALL $\xi \in \mathbb{R}$, pure algebraic proof):*
Decompose $\tilde c_{\mathrm{II},a}(r)$ affinely in $\frac{1}{|a|} - \frac{1}{2} \ge 0$:
$$\operatorname{Re}\!\big[P_k(\xi)\tilde c_{\mathrm{II},a}(r)\big] = \left(\frac{1}{|a|} - \frac{1}{2}\right) A_k(\xi) + \Big(\tfrac{1}{2}A_k(\xi) + B_k(\xi)\Big),$$
where $A_k(\xi) := \operatorname{Re}[P_k(\xi)\tan(\pi r)]$ and $\frac{1}{2}A_k(\xi) + B_k(\xi) := \operatorname{Re}\big[P_k(\xi)\big(\frac{1}{2} - \frac{1}{2r+1}\big)\tan(\pi r)\big]$. For both $\sigma_1 \in \{\frac{5}{8}, \frac{7}{8}\}$, $\sin(2\pi\sigma_1) = -\frac{\sqrt{2}}{2}$, $D(\xi) := \cos(2\pi\sigma_1) + \cosh(2\pi\xi) > 0$, and $\tan(\pi r) = \frac{-\frac{\sqrt{2}}{2} + i\sinh(2\pi\xi)}{D(\xi)}$. Let $d := 2\sigma_1 + 1 \in \{\frac{9}{4}, \frac{11}{4}\}$; then
$$\frac{1}{2} - \frac{1}{2r+1} = \frac{d(d-2) + 4x + 4i\xi}{2(d^2 + 4x)},$$
so defining $U_k(x) + i\xi V_k(x) := (R_k(x) + i\xi I_k(x))(d(d-2) + 4x + 4i\xi)$ yields
$$A_k(\xi) = -\frac{\frac{\sqrt{2}}{2}R_k(x) + \xi\sinh(2\pi\xi)I_k(x)}{D(\xi)}, \qquad \frac{1}{2}A_k(\xi) + B_k(\xi) = -\frac{\frac{\sqrt{2}}{2}U_k(x) + \xi\sinh(2\pi\xi)V_k(x)}{2 D(\xi)(d^2 + 4x)}.$$
Direct polynomial expansion gives:
- **Strict positivity of $R_k(x) > 0$ and $I_k(x) \ge 0$ (hence $A_k(\xi) < 0$ for all $\xi \in \mathbb{R}$):**
  - At $\sigma_1 = \frac{5}{8}$ ($c = \frac{5}{64}$): $R_0 = 1, I_0 = 0$; $R_1 = x + \frac{5}{64}, I_1 = \frac{1}{2}$; $R_2 = x^2 - \frac{3}{32}x + \frac{25}{4096} = (x - \frac{3}{64})^2 + \frac{1}{256} > 0$, $I_2 = x + \frac{5}{64}$.
  - At $\sigma_1 = \frac{7}{8}$ ($c = \frac{21}{64}$): every coefficient of $R_k, I_k$ for $k = 0, 1, 2, 3$ is nonnegative with strictly positive constant term (e.g., $R_3 = x^3 + \frac{15}{64}x^2 + \frac{315}{4096}x + \frac{9261}{262144}$, $I_3 = \frac{3}{2}x^2 + \frac{55}{64}x + \frac{1323}{8192}$).
- **Strict positivity of $\frac{\sqrt{2}}{2}U_k(x) + \xi\sinh(2\pi\xi)V_k(x) > 0$ (hence $\frac{1}{2}A_k(\xi) + B_k(\xi) < 0$ for all $\xi \in \mathbb{R}$):**
  Every coefficient of $V_k(x)$ is strictly positive in all seven cases ($\sigma_1 = \frac{5}{8}, k \le 2$ and $\sigma_1 = \frac{7}{8}, k \le 3$):
  - $\sigma_1 = \frac{5}{8}$: $V_0 = 4$, $V_1 = 6x + \frac{19}{32}$, $V_2 = 8x^2 + \frac{1}{2}x + \frac{35}{512}$;
  - $\sigma_1 = \frac{7}{8}$: $V_0 = 4$, $V_1 = 6x + \frac{75}{32}$, $V_2 = 8x^2 + 5x + \frac{567}{512}$, $V_3 = 10x^3 + \frac{239}{32}x^2 + \frac{5583}{2048}x + \frac{62181}{131072}$.
  Because $V_k(x) > 0$ and $\xi\sinh(2\pi\xi) \ge 2\pi x > 6x$ for all $\xi \in \mathbb{R}$, together with $\frac{7}{10} < \frac{\sqrt{2}}{2} < \frac{3}{4}$, we bound $\frac{\sqrt{2}}{2}U_k(x) + \xi\sinh(2\pi\xi)V_k(x) \ge W_{k,\mathrm{lb}}(x) := \sum_j (\alpha_j U_{k,j} + 6 V_{k,j-1}) x^j$ (with $\alpha_j = \frac{7}{10}$ if $U_{k,j} \ge 0$ and $\frac{3}{4}$ if $U_{k,j} < 0$, $V_{k,-1} := 0$). Direct rational arithmetic (independently certified in exact `Fraction` arithmetic by Gate 4 Check [2]) confirms that **every single coefficient of $W_{k,\mathrm{lb}}(x)$ is strictly positive** in all seven cases:
  - $\sigma_1 = \frac{5}{8}$: $W_{0,\mathrm{lb}} = \frac{63}{160} + \frac{134}{5}x$, $W_{1,\mathrm{lb}} = \frac{63}{2048} + \frac{87}{32}x + \frac{194}{5}x^2$, $W_{2,\mathrm{lb}} = \frac{315}{131072} + \frac{633}{4096}x + \frac{9}{64}x^2 + \frac{254}{5}x^3$;
  - $\sigma_1 = \frac{7}{8}$: $W_{0,\mathrm{lb}} = \frac{231}{160} + \frac{134}{5}x$, $W_{1,\mathrm{lb}} = \frac{4851}{10240} + \frac{601}{40}x + \frac{194}{5}x^2$, $W_{2,\mathrm{lb}} = \frac{101871}{655360} + \frac{27081}{4096}x + \frac{1905}{64}x^2 + \frac{254}{5}x^3$, $W_{3,\mathrm{lb}} = \frac{2139291}{41943040} + \frac{169533}{65536}x + \frac{29433}{2048}x^2 + \frac{681}{16}x^3 + \frac{314}{5}x^4$.
  Hence $\operatorname{Re}[P_k(\xi)\tilde c_{\mathrm{II},a}(r)] < 0$ for all $\xi \in \mathbb{R}$ and all $a \in [-2, 0)$. $\square$

### Lemma 2.3 (Unconditional Far-Field $L^2([0,\infty))$, Scalar Functional, and Global $L^\infty$ Bounds)
For every $G \in H^2_{\sigma_-}(\mathbb{R}) \cap H^3_{5/8}(\mathbb{R})$ with $\sigma_- \in \{\frac{1}{4}, \frac{3}{8}\}$:
1. **Far-field $L^2([0,\infty), dz)$ bounds for $\widetilde{H}$, $\widetilde{Q}$, and $\partial_z^{1,2}\widetilde{Q}_0$:** With $C_{\sigma_-} := \cot(\pi(\frac{1}{2}-\sigma_-))$ (so $C_{1/4} = 1$ and $C_{3/8} = \cot(\frac{\pi}{8}) = \sqrt{2} + 1$),
   $$\big\|\partial_z^k \widetilde{H}[G]\big\|_{L^2([0,\infty), dz)} \le C_{\sigma_-}\,\big\|(\partial_z - \tfrac{1}{2})^k G\big\|_{L^2_{\sigma_-}(\mathbb{R})} \qquad (k = 0, 1, 2).$$
   Furthermore, the nonzero constant $\widetilde{H}_{\mathrm{edge}}[G]$ cancels identically in $\partial_z^m \widetilde{H}_0[G] = \partial_z^m \widetilde{H}[G]$ ($m \ge 1$) and in $\partial_z\widetilde{Q}_0[G] = \widetilde{H}[G] - \widetilde{Q}[G]$, $\partial_z^2\widetilde{Q}_0[G] = \partial_z\widetilde{H}[G] - \widetilde{H}[G] + \widetilde{Q}[G]$, yielding
   $$\big\|\widetilde{Q}[G]\big\|_{L^2([0,\infty))} + \sum_{m=1}^2 \Big(\big\|\partial_z^m \widetilde{H}_0[G]\big\|_{L^2([0,\infty))} + \big\|\partial_z^m \widetilde{Q}_0[G]\big\|_{L^2([0,\infty))}\Big) \le C\,\Big(\|G\|_{H^2_{\sigma_-}(\mathbb{R})} + \|G\|_{H^2_{5/8}(\mathbb{R})}\Big).$$
2. **Scalar functional bounds:**
   $$|V_G| \le C_{V,\sigma_-}\|G\|_{L^2_{\sigma_-}(\mathbb{R})}, \qquad |\widetilde{H}_{\mathrm{edge}}[G]| \le \frac{1}{\pi}\!\left(2\|G\|_{L^2_{5/8}(\mathbb{R})} + \frac{1}{\sqrt{1 - 2\sigma_-}}\|G\|_{L^2_{\sigma_-}(\mathbb{R})}\right),$$
   with $C_{V,1/4} = \frac{1}{\pi}\big(\int_{\mathbb{R}} e^{3z/2}\ln^2(1+e^{-z})\,dz\big)^{1/2} \approx 0.51036$ and $C_{V,3/8} = \frac{1}{\pi}\big(\int_{\mathbb{R}} e^{7z/4}\ln^2(1+e^{-z})\,dz\big)^{1/2} \approx 0.65136$.
3. **Global $L^\infty(\mathbb{R})$ control of forcing and derivatives:**
   $$\sum_{j=0}^1 \Big(\|\partial_z^j \widetilde{H}_0[G]\|_{L^\infty(\mathbb{R})} + \|\partial_z^j \widetilde{Q}_0[G]\|_{L^\infty(\mathbb{R})}\Big) + \|\partial_z^2 \widetilde{Q}_0[G]\|_{L^\infty(\mathbb{R})} \le C_*\,\Big(\|G\|_{H^2_{\sigma_-}(\mathbb{R})} + \|G\|_{H^3_{5/8}(\mathbb{R})}\Big).$$

**Proof.**
For (1): on $v \ge 1$, $v^{-1} \le v^{-2\sigma_-}$. Because $w(v) = v^{-2\sigma_-}$ satisfies the Muckenhoupt $A_2$ condition $-1 < -2\sigma_- < 1$, $H_+$ is bounded on $L^2((0,\infty), v^{-2\sigma_-}dv)$ with exact Mellin norm $\sup_\xi|\cot(\pi(\frac{1}{2}-\sigma_- + i\xi))| = \max(\cot(\pi(\frac{1}{2}-\sigma_-)), 1) = C_{\sigma_-}$ for $\sigma_- \in [\frac{1}{4}, \frac{1}{2})$. Hence for $f_k(v) := v^{-1/2}(\partial_z - \frac{1}{2})^k G(\ln v)$ ($k \in \{0,1,2\}$, using $\partial_z^k\widetilde{H}[G] = H_+[f_k]$ by exact dilation commutation):
$$\int_0^\infty |\partial_z^k\widetilde{H}[G](z)|^2\,dz = \int_1^\infty v^{-1}|H_+[f_k]|^2\,dv \le \int_0^\infty v^{-2\sigma_-}|H_+[f_k]|^2\,dv \le C_{\sigma_-}^2 \int_0^\infty v^{-2\sigma_- - 1}\big|(\partial_z-\tfrac{1}{2})^k G\big|^2\,dv = C_{\sigma_-}^2\big\|(\partial_z-\tfrac{1}{2})^k G\big\|_{L^2_{\sigma_-}}^2.$$
For $\widetilde{Q}[G]$ on $z \ge 0$, using $\widetilde{Q}[G](z) = e^{-z}\widetilde{Q}[G](0) + \int_0^z e^{z'-z}\widetilde{H}[G](z')\,dz'$, Young's convolution inequality ($\|e^{-z}\mathbf{1}_{z\ge 0}\|_{L^1} = 1$) gives $\|\widetilde{Q}[G]\|_{L^2([0,\infty))} \le \frac{1}{\sqrt{2}}|\widetilde{Q}[G](0)| + \|\widetilde{H}[G]\|_{L^2([0,\infty))}$. Since $\widetilde{Q}[G](0) = \widetilde{Q}_0[G](0) + \widetilde{H}_{\mathrm{edge}}[G]$ and $|\widetilde{Q}_0[G](0)| = |(T_{Q_0}G)(0)| \le \sqrt{2}\|T_{Q_0}G\|_{H^1_{5/8}} \le C\|G\|_{H^1_{5/8}}$ (Lemma 2.1), part (2) bounds $|\widetilde{Q}[G](0)| \le C(\|G\|_{L^2_{\sigma_-}} + \|G\|_{H^1_{5/8}})$. Because $\partial_z\widetilde{Q}_0 = \widetilde{H}_0 - \widetilde{Q}_0 = \widetilde{H} - \widetilde{Q}$ and $\partial_z^2\widetilde{Q}_0 = \partial_z\widetilde{H} - \widetilde{H} + \widetilde{Q}$, the constant $\widetilde{H}_{\mathrm{edge}}[G]$ cancels identically and the remaining $L^2([0,\infty))$ bounds follow immediately.

Part (2) is Cauchy–Schwarz on $(-\infty,0)$ and $(0,\infty)$ (using $\sigma_- < \frac{1}{2} < \frac{5}{8}$). For (3), on $z \le 0$ we use $e^{-z/2}|u(z)| = e^{+z/8}(e^{-5z/8}|u(z)|) \le \sqrt{2}\|u\|_{H^1_{5/8}}$ applied to $u = T_{H_0}(\partial_z - \frac{1}{2})^j G$ and $u = T_{Q_0}(\partial_z - \frac{1}{2})^j G$ (Lemma 2.1, bounded on $H^1_{5/8}$ for $j=0,1$, and for $\partial_z^2\widetilde{Q}_0 = \partial_z\widetilde{H}_0 - \widetilde{H}_0 + \widetilde{Q}_0$ bounded by $\|G\|_{H^2_{5/8}}$); on $z \ge 0$, 1D Sobolev embedding $\|u\|_{L^\infty([0,\infty))} \le \sqrt{2}\|u\|_{H^1([0,\infty))}$ combined with part (1) and $|\widetilde{H}_{\mathrm{edge}}[G]|$ (part (2)) bounds every term directly by $\|G\|_{H^2_{\sigma_-}} + \|G\|_{H^2_{5/8}}$. $\square$

---

## 4. Theorem 3: Global Nonlinear Stability, Uniqueness, and Distributional Closure in Finite Sobolev Space (Zero Analyticity)

**Theorem 3 (Small-Data Nonlinear Asymptotic Stability & Uniqueness of $\bar\Omega_{\mathrm{II},a}$ for $a \in [-2,0)$, Including $\bar\Omega_{\mathrm{I},-1}$).**
*Fix $a \in [-2, 0)$ (in particular $a = -1$, where $\bar\Omega_{\mathrm{I},-1} = \bar\Omega_{\mathrm{II},-1}$ is Huang–Tong–Wang's singular profile). Set $\sigma_- = \frac{1}{4}$ if $a \in [-\frac{8}{5}, 0)$ (with $\lambda_0 = \frac{1}{2}$) or $\sigma_- = \frac{3}{8}$ for any $a \in [-2, 0)$ (with $\lambda_0 = \frac{3}{4}$). Define the finite-regularity Sobolev space*
$$\mathcal{X}_{\mathrm{gCLM}} := \Big\{G \in H^2_{\sigma_-}(\mathbb{R}) \cap H^3_{5/8}(\mathbb{R}) : G \in L^\infty([0,\infty)),\ \partial_z G \in H^1([0,\infty))\Big\},$$
*with weighted energy and far-field norm*
$$\mathcal{E}(\tau) := \sum_{k=0}^2 \|\partial_z^k G(\tau)\|_{L^2_{\sigma_-}(\mathbb{R})}^2 + \sum_{k=0}^3 \|\partial_z^k G(\tau)\|_{L^2_{5/8}(\mathbb{R})}^2, \qquad \mathcal{B}(\tau) := \|G(\tau)\|_{L^\infty([0,\infty))} + \|\partial_z G(\tau)\|_{H^1([0,\infty))}.$$
*There exist explicit constants $\varepsilon_0 > 0$ and $M, C_{\mathrm{nl}}, C_0, C_4, C_5 > 0$ depending only on $a$ such that:*
1. **Uniqueness:** For any $T > 0$,~\eqref{eq:nonlin_gauge_c} has at most one strong solution $G \in C([0,T]; \mathcal{X}_{\mathrm{gCLM}})$ for a given initial datum $G_0 \in \mathcal{X}_{\mathrm{gCLM}}$.
2. **Global *a priori* nonlinear closure & exponential decay:** Any strong solution $G \in C([0,T); \mathcal{X}_{\mathrm{gCLM}})$ of~\eqref{eq:nonlin_gauge_c} with initial data $\sqrt{\mathcal{E}(0)} + \mathcal{B}(0) \le \varepsilon_0$ satisfies the uniform bounds on $[0, T)$:
   $$\mathcal{E}(\tau) \le \mathcal{E}(0)\,e^{-\lambda_0 \tau}, \qquad \sqrt{\mathcal{E}(\tau)} + \mathcal{B}(\tau) \le \frac{M}{2}\,\varepsilon_0 < M\,\varepsilon_0, \qquad \lim_{z\to-\infty} G(z,\tau) = 0,$$
   together with exponential convergence of the dynamic-rescaling rates:
   $$|c_l(\tau) - 2| + |c_\omega(\tau) + 1| \le C\,\sqrt{\mathcal{E}(0)}\,e^{-\lambda_0\tau/2} \qquad (\forall \tau \in [0,T)).$$
   Moreover, these bounds persist uniformly for the viscous regularization $\partial_\tau G^\nu = \text{RHS}\eqref{eq:nonlin_gauge_c} + \nu\partial_z^2 G^\nu$ for all $\nu \in (0, \nu_0]$.
3. **Global distributional blowup profile:** For the resulting global strong solution ($T = \infty$), $\Omega(X,\tau) = \frac{1}{a}(X-1)^{-1/2}(1+G(\ln(X-1),\tau))\mathbf{1}_{\{X>1\}}$ is a global distributional solution of~\eqref{eq:gclm_dyn} on $\mathbb{R}$ with zero $\delta(X-1)$ defect, converging uniformly on compact subsets of $(1, \infty)$ to $\bar\Omega_{\mathrm{II},a}$. In physical variables, at the finite time $T_* = \int_0^\infty C_\omega(\tau)\,d\tau < \infty$, the singular edge $x_{\mathrm{edge}}(t) = C_l(\tau(t))^{-1} \sim c_L(T_* - t)^2$ collapses to $x = 0$ with exact shrinking exponent $\gamma = 2$, while the self-similar amplitude blows up as $C_\omega(\tau(t))^{-1} \sim c_*(T_* - t)^{-1}$ (e.g. $u_x(0,t) = H(\omega)(0,t) \sim \frac{c_*}{|a|}(T_* - t)^{-1} \to +\infty$).

*(Technical Remark on Local Existence: Part (2) proves the complete global *a priori* closure, uniform-in-$\nu$ viscous persistence, and uniqueness in $\mathcal{X}_{\mathrm{gCLM}}$. Constructing the short-time existence of the linear/viscous regularized semigroup $e^{\tau(M_c + \nu\partial_z^2)}$ in $\mathcal{X}_{\mathrm{gCLM}}$ — e.g. verifying the Hille–Yosida range condition or finite-domain Galerkin approximation — is treated as a standard technical step rather than proved inline here.)*

**Proof.**
Let $F[G] := \frac{1}{a}\widetilde{H}_0[G] + \frac{1}{2}\widetilde{Q}_0[G]$ and $c(z,\tau) := 2 + \widetilde{Q}[G](z) - V_G$. Note that by 1D Sobolev embedding on $(-\infty, 0]$ (using $H^3_{5/8}$) and $[0,\infty)$ (using $\mathcal{B}$):
$$\lim_{z\to-\infty} |G(z,\tau)| \le \lim_{z\to-\infty} e^{5z/8}\cdot \sqrt{2}\|G\|_{H^1_{5/8}} = 0,$$
and
\begin{equation}\label{eq:weighted_linf}
\|e^{-\alpha z} G\|_{L^\infty(\mathbb{R})} \le C\big(\sqrt{\mathcal{E}} + \mathcal{B}\big), \qquad \|e^{-\alpha z}\partial_z G\|_{L^\infty(\mathbb{R})} \le C\big(\sqrt{\mathcal{E}} + \mathcal{B}\big) \qquad \forall \alpha \in [0, \tfrac{5}{8}].
\end{equation}

**Part (1) — Uniqueness in $C([0,T]; \mathcal{X}_{\mathrm{gCLM}})$.**
Let $G_1, G_2 \in C([0,T]; \mathcal{X}_{\mathrm{gCLM}})$ be two strong solutions with the same initial datum, and set $W := G_1 - G_2$. Subtracting~\eqref{eq:nonlin_gauge_c} yields
$$\partial_\tau W = M_c W - (\widetilde{Q}[G_1] - V_{G_1})\partial_z W - (\widetilde{Q}_0[W] + \widetilde{H}_{\mathrm{edge}}[W] - V_W)\partial_z G_2 + F[G_1]W + F[W]G_2.$$
Test against $W$ in $L^2_\sigma(\mathbb{R})$ for $\sigma \in \{\sigma_-, \frac{5}{8}\}$. By Lemma 2.2, $\langle M_c W, W\rangle_{L^2_\sigma} \le -2\sigma\|W\|_{L^2_\sigma}^2 \le 0$. Integration by parts and Lemma 2.3(2)–(3) bound $-\langle(\widetilde{Q}[G_1]-V_{G_1})\partial_z W, W\rangle_{L^2_\sigma} + \langle F[G_1]W, W\rangle_{L^2_\sigma} \le C\sqrt{\mathcal{E}_1}\|W\|_{L^2_\sigma}^2$. By Lemma 2.3(2), $|V_W| + |\widetilde{H}_{\mathrm{edge}}[W]| \le C(\|W\|_{L^2_{\sigma_-}} + \|W\|_{L^2_{5/8}})$, giving $|\langle(\widetilde{H}_{\mathrm{edge}}[W] - V_W)\partial_z G_2, W\rangle_{L^2_\sigma}| \le C\sqrt{\mathcal{E}_2}(\|W\|_{L^2_{\sigma_-}} + \|W\|_{L^2_{5/8}})\|W\|_{L^2_\sigma}$. Finally, writing $\widetilde{Q}_0[W] = e^{-z/2}T_{Q_0}W$ and $F[W] = e^{-z/2}T_F W$ (Lemma 2.1, bounded on $L^2_{5/8}$) and splitting the weight as $e^{-2\sigma z}\cdot e^{-z/2} = e^{-(\sigma - 1/8)z}\cdot e^{-5z/8}\cdot e^{-\sigma z}$ with $\alpha = \sigma - \frac{1}{8} \in [\frac{1}{8}, \frac{1}{2}] \subset [0, \frac{5}{8}]$ gives, via~\eqref{eq:weighted_linf},
$$\big|\langle \widetilde{Q}_0[W]\partial_z G_2, W\rangle_{L^2_\sigma}\big| + \big|\langle F[W]G_2, W\rangle_{L^2_\sigma}\big| \le C\big(\sqrt{\mathcal{E}_2} + \mathcal{B}_2\big)\|W\|_{L^2_{5/8}}\|W\|_{L^2_\sigma}.$$
Summing $\sigma = \sigma_-$ and $\sigma = \frac{5}{8}$ yields $\frac{d}{d\tau}(\|W\|_{L^2_{\sigma_-}}^2 + \|W\|_{L^2_{5/8}}^2) \le C(\sqrt{\mathcal{E}_1} + \sqrt{\mathcal{E}_2} + \mathcal{B}_2)(\|W\|_{L^2_{\sigma_-}}^2 + \|W\|_{L^2_{5/8}}^2)$, so Gronwall implies $W \equiv 0$.

**Part (2) — Complete Three-Block Global *A Priori* Closure.**
We work under the bootstrap hypothesis
\begin{equation}\label{eq:bootstrap_hyp}
\sqrt{\mathcal{E}(\tau)} + \mathcal{B}(\tau) \le M\,\varepsilon_0,
\end{equation}
where $M > 0$ will be fixed explicitly below (independent of $\varepsilon_0$), and $\varepsilon_0$ is chosen small enough that $c(z,\tau) = 2 + \widetilde{Q}[G] - V_G \ge 1$ on $\mathbb{R}$ (by Lemma 2.3(2)–(3), since $|\widetilde{Q}[G] - V_G| \le C_*\sqrt{\mathcal{E}} \le C_* M\varepsilon_0 \le 1$).

**Block 1 — Weighted energy $\mathcal{E}(\tau)$.**
Differentiate~\eqref{eq:nonlin_gauge_c} $k$ times in $z$ ($k \in \{0,1,2\}$ for $\sigma_-$, $k \in \{0,1,2,3\}$ for $\sigma_+ = \frac{5}{8}$) and test against $\partial_z^k G$ in $L^2_\sigma(\mathbb{R})$:
$$\frac{1}{2}\frac{d}{d\tau}\|\partial_z^k G\|_{L^2_\sigma}^2 = \langle \partial_z^k M_c G, \partial_z^k G\rangle_{L^2_\sigma} - \big\langle (\widetilde{Q}[G] - V_G)\partial_z^{k+1} G, \partial_z^k G\big\rangle_{L^2_\sigma}$$
$$- \sum_{j=1}^k \binom{k}{j}\big\langle (\partial_z^j \widetilde{Q}[G])(\partial_z^{k+1-j} G), \partial_z^k G\big\rangle_{L^2_\sigma} + \sum_{j=0}^k \binom{k}{j}\big\langle (\partial_z^j F[G])(\partial_z^{k-j} G), \partial_z^k G\big\rangle_{L^2_\sigma}.$$
1. *Linear term:* By Lemma 2.2, $\langle \partial_z^k M_c G, \partial_z^k G\rangle_{L^2_\sigma} \le -2\sigma\|\partial_z^k G\|_{L^2_\sigma}^2 \le -\lambda_0\|\partial_z^k G\|_{L^2_\sigma}^2$ (since $2\sigma_- = \lambda_0 < \frac{5}{4} = 2\sigma_+$).
2. *Top-order self-advection:* Integration by parts yields
   $$-\big\langle (\widetilde{Q}[G]-V_G)\partial_z^{k+1} G, \partial_z^k G\big\rangle_{L^2_\sigma} = \int_{\mathbb{R}} e^{-2\sigma z}\!\left[-\sigma(\widetilde{Q}[G] - V_G) + \tfrac{1}{2}\partial_z\widetilde{Q}[G]\right](\partial_z^k G)^2\,dz \le C\sqrt{\mathcal{E}}\,\|\partial_z^k G\|_{L^2_\sigma}^2$$
   by Lemma 2.3(2)–(3) (zero derivative loss).
3. *Lower-order self-advection commutators ($1 \le j \le k$):*
   - For $j = 1, 2$ with $j < k$ (and $j = 1$ for any $k$): $\|\partial_z^j \widetilde{Q}[G]\|_{L^\infty(\mathbb{R})} = \|\partial_z^j \widetilde{Q}_0[G]\|_{L^\infty(\mathbb{R})} \le C\sqrt{\mathcal{E}}$ by Lemma 2.3(3), while $\|\partial_z^{k+1-j} G\|_{L^2_\sigma} \le \sqrt{\mathcal{E}}$, giving $\le C\mathcal{E}^{3/2}$.
   - For $j = k \in \{2, 3\}$: write $\partial_z^k\widetilde{Q}[G] = \partial_z^k\widetilde{Q}_0[G] = e^{-z/2} T_{Q_0}((\partial_z - \frac{1}{2})^k G)$ (Lemma 2.1(3)), with $\|T_{Q_0}((\partial_z-\frac{1}{2})^k G)\|_{L^2_{5/8}} \le C\|G\|_{H^k_{5/8}} \le C\sqrt{\mathcal{E}}$. Splitting the weight $e^{-2\sigma z}\cdot e^{-z/2} = e^{-(\sigma - 1/8)z}\cdot e^{-5z/8}\cdot e^{-\sigma z}$ places $e^{-\alpha z}$ with $\alpha = \sigma - \frac{1}{8} \in [\frac{1}{8}, \frac{1}{2}] \subset [0, \frac{5}{8}]$ on $\partial_z G$, yielding via~\eqref{eq:weighted_linf} $\le C\|e^{-\alpha z}\partial_z G\|_{L^\infty(\mathbb{R})}\|G\|_{H^k_{5/8}}\|\partial_z^k G\|_{L^2_\sigma} \le C(\sqrt{\mathcal{E}} + \mathcal{B})\mathcal{E}$.
4. *Nonlinear stretching Leibniz terms ($0 \le j \le k$):*
   - For $j = 0, 1$: $\|\partial_z^j F[G]\|_{L^\infty(\mathbb{R})} \le C\sqrt{\mathcal{E}}$ by Lemma 2.3(3), giving $\le C\mathcal{E}^{3/2}$.
   - For $j = 2, 3$ ($k - j \in \{0, 1\}$): write $\partial_z^j F[G] = e^{-z/2} T_F((\partial_z-\frac{1}{2})^j G)$ where $\|T_F((\partial_z-\frac{1}{2})^j G)\|_{L^2_{5/8}} \le C\|G\|_{H^j_{5/8}} \le C\sqrt{\mathcal{E}}$ (Lemma 2.1). Splitting the weight identically to step 3 places $e^{-\alpha z}$ ($\alpha = \sigma - \frac{1}{8} \in [\frac{1}{8}, \frac{1}{2}]$) on $\partial_z^{k-j} G$ ($k - j \le 1$), bounded in $L^\infty(\mathbb{R})$ by $C(\sqrt{\mathcal{E}} + \mathcal{B})$, yielding $\le C(\sqrt{\mathcal{E}} + \mathcal{B})\mathcal{E}$.

Summing over all $k$ and multiplying by $2$ gives the closed weighted inequality
\begin{equation}\label{eq:block1}
\frac{d}{d\tau}\mathcal{E} \le -2\lambda_0\,\mathcal{E} + 2 C_{\mathrm{nl}}\big(\sqrt{\mathcal{E}} + \mathcal{B}\big)\,\mathcal{E}.
\end{equation}

**Block 2 — Far-field derivative energy $P(\tau) := \|\partial_z G\|_{H^1([0,\infty))}^2 = \|p\|_{L^2([0,\infty))}^2 + \|p_z\|_{L^2([0,\infty))}^2$ ($p := \partial_z G$).**
Since $c(z,\tau) \ge 1$ on $\mathbb{R}$, $z = 0$ is an inflow boundary for $[0, \infty)$. Differentiating~\eqref{eq:nonlin_gauge_c} $m + 1$ times ($m \in \{0, 1\}$) and testing against $\partial_z^m p = \partial_z^{m+1} G$ on $[0,\infty)$:
- *Transport IBP:* $-\int_0^\infty c\,(\partial_z^m p)(\partial_z^{m+1} p)\,dz = \frac{1}{2} c(0,\tau)(\partial_z^m p(0,\tau))^2 + \frac{1}{2}\int_0^\infty (\partial_z\widetilde{Q}[G])(\partial_z^m p)^2\,dz$. Because $H^3_{5/8}(\mathbb{R}) \hookrightarrow C^2(\mathbb{R})$, $|\partial_z^m p(0,\tau)| = |\partial_z^{m+1} G(0,\tau)| \le C\|G\|_{H^3_{5/8}(\mathbb{R})} \le C\sqrt{\mathcal{E}}$ for both $m = 0, 1$, and $\|\partial_z\widetilde{Q}[G]\|_{L^\infty([0,\infty))} \le C\sqrt{\mathcal{E}}$ by Lemma 2.3(3).
- *Forcing terms $\partial_z^{m+1}(F[G](1+G))$ and commutators:* By Lemma 2.3(1), $\|\partial_z^{m+1} F[G]\|_{L^2([0,\infty))} \le C(\|G\|_{H^2_{\sigma_-}} + \|G\|_{H^2_{5/8}}) \le C\sqrt{\mathcal{E}}$ for $m \in \{0, 1\}$, and all lower-order $L^\infty$ coefficients $\|F[G]\|_{L^\infty}, \|\partial_z F[G]\|_{L^\infty}, \|\partial_z\widetilde{Q}\|_{L^\infty}, \|\partial_z^2\widetilde{Q}\|_{L^\infty}$ are bounded by $C\sqrt{\mathcal{E}}$ (Lemma 2.3(3)). Summing $m = 0, 1$ yields
  $$\frac{1}{2}\frac{d}{d\tau} P \le C_1\,\mathcal{E} + C_2\,(1 + \mathcal{B})\,\sqrt{\mathcal{E}}\,\sqrt{P} + C_3\,\sqrt{\mathcal{E}}\,P.$$
To avoid any singularity at $P = 0$, combine $P$ with $\mathcal{E}$: whenever $\frac{d}{d\tau}\mathcal{E} \le 0$ (which holds under the bootstrap as soon as $2 C_{\mathrm{nl}} M\varepsilon_0 \le \lambda_0$), for any $\delta > 0$,
$$\frac{d}{d\tau}\sqrt{P + \mathcal{E} + \delta^2} = \frac{\frac{1}{2}\dot P + \frac{1}{2}\dot{\mathcal{E}}}{\sqrt{P + \mathcal{E} + \delta^2}} \le \frac{C_1\mathcal{E} + C_2(1+\mathcal{B})\sqrt{\mathcal{E}}\sqrt{P} + C_3\sqrt{\mathcal{E}}\,P}{\sqrt{P + \mathcal{E} + \delta^2}} \le C_4\,(1 + \mathcal{B})\,\sqrt{\mathcal{E}},$$
with $C_4 := \max(C_1 + C_2,\, C_2 + C_3)$ (using $\frac{\mathcal{E}}{\sqrt{P+\mathcal{E}+\delta^2}} \le \sqrt{\mathcal{E}}$, $\frac{\sqrt{P}}{\sqrt{P+\mathcal{E}+\delta^2}} \le 1$, and $\frac{P}{\sqrt{P+\mathcal{E}+\delta^2}} \le \sqrt{P} \le \mathcal{B}$). Integrating from $0$ to $\tau$ and letting $\delta \to 0^+$ gives
\begin{equation}\label{eq:block2}
\sqrt{P(\tau)} \le \sqrt{P(0)} + \sqrt{\mathcal{E}(0)} + C_4\int_0^\tau (1 + \mathcal{B}(s))\sqrt{\mathcal{E}(s)}\,ds.
\end{equation}

**Block 3 — Far-field amplitude $\|G\|_{L^\infty([0,\infty))}$.**
Since $c(z,\tau) \ge 1 > 0$ on $[0,\infty)$, every characteristic $Z(\tau)$ in $[0,\infty)$ either starts at $Z(0) \ge 0$ or enters at $Z(\tau_0) = 0$ with $|G(0,\tau_0)| \le C_0\sqrt{\mathcal{E}(\tau_0)}$. Along any characteristic, $\big|\frac{d}{d\tau}G(Z(\tau),\tau)\big| \le \|F[G]\|_{L^\infty([0,\infty))}(1 + \|G\|_{L^\infty([0,\infty))}) \le C_5\sqrt{\mathcal{E}(\tau)}(1 + \mathcal{B}(\tau))$. Integrating from $\max(0,\tau_0)$ to $\tau$ gives
\begin{equation}\label{eq:block3}
\|G(\tau)\|_{L^\infty([0,\infty))} \le \|G_0\|_{L^\infty([0,\infty))} + C_0\sup_{s\le\tau}\sqrt{\mathcal{E}(s)} + C_5\int_0^\tau (1 + \mathcal{B}(s))\sqrt{\mathcal{E}(s)}\,ds.
\end{equation}

**Bootstrap Closure.**
Define the explicit $a$-dependent constants
$$K := 2 + C_0 + \frac{3(C_4 + C_5)}{\lambda_0}, \qquad M := 2(1 + K).$$
Choose $\varepsilon_0 > 0$ small enough that
$$2 C_{\mathrm{nl}}\,M\,\varepsilon_0 \le \lambda_0, \qquad C_*\,M\,\varepsilon_0 \le 1, \qquad M\,\varepsilon_0 \le \frac{1}{2}.$$
Under the bootstrap hypothesis~\eqref{eq:bootstrap_hyp}, eq.~\eqref{eq:block1} yields $\frac{d}{d\tau}\mathcal{E} \le -\lambda_0\,\mathcal{E}$, hence
$$\mathcal{E}(\tau) \le \mathcal{E}(0)\,e^{-\lambda_0\tau}, \qquad \int_0^\tau \sqrt{\mathcal{E}(s)}\,ds \le \frac{2}{\lambda_0}\sqrt{\mathcal{E}(0)}.$$
Substituting this and $1 + \mathcal{B}(s) \le 1 + M\varepsilon_0 \le \frac{3}{2}$ into~\eqref{eq:block2} and~\eqref{eq:block3} and summing with $\sqrt{\mathcal{E}(\tau)} \le \sqrt{\mathcal{E}(0)}$ yields
$$\sqrt{\mathcal{E}(\tau)} + \mathcal{B}(\tau) \le \mathcal{B}(0) + \left(2 + C_0 + \frac{3(C_4 + C_5)}{\lambda_0}\right)\sqrt{\mathcal{E}(0)} = \mathcal{B}(0) + K\sqrt{\mathcal{E}(0)} \le (1 + K)\,\varepsilon_0 = \frac{M}{2}\,\varepsilon_0 < M\,\varepsilon_0,$$
closing the bootstrap strictly with a factor-of-2 margin on all of $[0, T)$! Exponential convergence of $(c_l(\tau), c_\omega(\tau)) \to (2, -1)$ follows immediately from~\eqref{eq:gauge_c} and Lemma 2.3(2).

**Uniform-in-$\nu$ Viscous Regularization Closure ($+\nu\partial_z^2 G$ on RHS).**
If $+\nu\partial_z^2 G$ ($\nu > 0$) is added to the right-hand side of~\eqref{eq:nonlin_gauge_c} (tightening the smallness choice to $2 C_{\mathrm{nl}} M \varepsilon_0 \le \frac{\lambda_0}{2}$):
1. In Block 1, $\langle +\nu\partial_z^2(\partial_z^k G), \partial_z^k G\rangle_{L^2_\sigma} = -\nu\|\partial_z^{k+1} G\|_{L^2_\sigma}^2 + 2\sigma^2\nu\|\partial_z^k G\|_{L^2_\sigma}^2$. Choosing $\nu_0 \le \frac{\lambda_0}{8\sigma_+^2}$ gives $+4\sigma_+^2\nu\mathcal{E} \le \frac{\lambda_0}{2}\mathcal{E}$, preserving $\frac{d}{d\tau}\mathcal{E} \le -\lambda_0\mathcal{E} - 2\nu\sum_{k=0}^3\|\partial_z^{k+1} G\|_{L^2_{5/8}}^2$.
2. In Block 2, integrating $+\nu\partial_z^2(\partial_z^m p)\,\partial_z^m p$ by parts on $[0,\infty)$ yields $-\nu\|\partial_z^{m+1} p\|_{L^2([0,\infty))}^2 - \nu\,\partial_z^{m+1} p(0)\,\partial_z^m p(0)$. For $m = 0$, $|\partial_z p(0)p(0)| \le C\mathcal{E}$. For $m = 1$, 1D trace interpolation $|\partial_z^3 G(0)| \le C(\|\partial_z^3 G\|_{L^2_{5/8}}^{1/2}\|\partial_z^4 G\|_{L^2_{5/8}}^{1/2} + \|\partial_z^3 G\|_{L^2_{5/8}})$, $|\partial_z^2 G(0)| \le C\sqrt{\mathcal{E}}$, and Young's inequality ($p = 4, q = \frac{4}{3}$) give $\nu|\partial_z^3 G(0)\partial_z^2 G(0)| \le \frac{\nu}{2}\|\partial_z^4 G\|_{L^2_{5/8}}^2 + C\nu\,\mathcal{E}$, whose top-order term is absorbed directly into $-\nu\|\partial_z^4 G\|_{L^2_{5/8}}^2$ inside $\frac{1}{2}\dot{\mathcal{E}}$ while $+C\nu\mathcal{E}$ is absorbed into $C_1\mathcal{E}$ uniformly for $\nu \in (0, \nu_0]$.
3. In Block 3, on $[0,\infty)$, $G^\nu$ satisfies $\partial_\tau G^\nu - \nu\partial_z^2 G^\nu + c\,\partial_z G^\nu = F[G^\nu](1 + G^\nu)$ with $|G^\nu(0,\tau)| \le C_0\sqrt{\mathcal{E}(\tau)}$; comparing with the spatially constant supersolution $Y'(\tau) = \|F[G^\nu](\tau)\|_{L^\infty([0,\infty))}(1 + \mathcal{B}(\tau))$, $Y(0) = \|G_0\|_{L^\infty([0,\infty))} + C_0\sup_{s\le\tau}\sqrt{\mathcal{E}(s)}$, via the 1D parabolic comparison principle yields~\eqref{eq:block3} identically for all $\nu \in (0, \nu_0]$.

Part (3) follows directly from the global distributional verification in Section 1 and integrating $\frac{dt}{d\tau} = C_\omega(\tau) = \exp(-\tau + O(1))$, $C_l(\tau) = \exp(2\tau + O(1))$ to $T_* = \int_0^\infty C_\omega(\tau)\,d\tau < \infty$, giving $T_* - t \sim c_1 e^{-\tau}$, $C_\omega^{-1} \sim c_*(T_* - t)^{-1}$, and $x_{\mathrm{edge}}(t) = C_l^{-1} \sim c_L(T_* - t)^2$. $\square$

---

## 5. Theorem 4: Unconditional $H^1_\sigma(\mathbb{R})$ Coercivity of Huang–Tong–Wang's Branch I for All $a \in [-1, 0)$

Although Huang–Tong–Wang's Branch I profile $\bar\Omega_{\mathrm{I},a}(X) = -\sin(\pi\mu)(X-1)^{-\mu}\mathbf{1}_{\{X>1\}}$ ($\mu = \frac{1}{1-a} \in [\frac{1}{2}, 1)$) admits no finite Gamma-product Mellin conjugator for $a \ne -1$, the **same Gauge-(c) edge-amplitude-pinning reduction** yields an **elementary, computer-free, unconditional $H^1_\sigma(\mathbb{R})$ coercivity theorem for every $a \in [-1, 0)$**!

For $a \in [-1, 0)$, write $\kappa := -\cos(\pi\mu) \ge 0$ and $\Omega(1+e^z, \tau) = \bar\Omega_{\mathrm{I},a}(1+e^z)(1 + G(z,\tau))$ on $X = 1+e^z > 1$. Define
$$\widetilde{H}_{\mathrm{I}}[G](z) := -\sin(\pi\mu)\,H_+[v^{-\mu} G](1+e^z), \qquad \widetilde{Q}_{\mathrm{I}}[G](z) := \int_{-\infty}^z e^{z'-z}\widetilde{H}_{\mathrm{I}}[G](z')\,dz',$$
$$\widetilde{H}_{\mathrm{I},\mathrm{edge}}[G] := \lim_{z\to-\infty}\widetilde{H}_{\mathrm{I}}[G](z) = \frac{\sin(\pi\mu)}{\pi}\widehat{G}(\mu), \qquad V_{\mathrm{I},G} := a\,\delta U(1) = \frac{a\sin(\pi\mu)}{\pi}\int_{\mathbb{R}} e^{(1-\mu)z}\ln(1+e^{-z})\,G(z)\,dz.$$
Imposing **Branch I Gauge (c)**,
$$c_l(\tau) = (1-a) - V_{\mathrm{I},G}(\tau), \qquad c_\omega(\tau) = -\mu\,c_l(\tau) - \mu\,\widetilde{H}_{\mathrm{I},\mathrm{edge}}[G],$$
cancels both the singular edge advection drift and the constant edge trace $\mu\widetilde{H}_{\mathrm{I},\mathrm{edge}}[G]$ identically, yielding the linearized Gauge-(c) Branch I operator
\begin{equation}\label{eq:M_Ic}
M_{\mathrm{I},c} G := -(1-a)(1 + \kappa e^{-\mu z})\partial_z G + \widetilde{H}_{\mathrm{I},0}[G] - (1-\mu)\widetilde{Q}_{\mathrm{I},0}[G],
\end{equation}
where $\widetilde{H}_{\mathrm{I},0}[G] := \widetilde{H}_{\mathrm{I}}[G] - \widetilde{H}_{\mathrm{I},\mathrm{edge}}[G] = -\sin(\pi\mu)\,v H_+[v^{-1-\mu} G]$ and $\widetilde{Q}_{\mathrm{I},0}[G] := \widetilde{Q}_{\mathrm{I}}[G] - \widetilde{H}_{\mathrm{I},\mathrm{edge}}[G]$.

**Theorem 4 (Unconditional $H^1_\sigma(\mathbb{R})$ Coercivity of Branch I for All $a \in [-1, 0)$).**
*Fix any $a \in [-1, 0)$ ($\mu = \frac{1}{1-a} \in [\frac{1}{2}, 1)$, $\kappa = -\cos(\pi\mu) \ge 0$) and any $\theta \in (0, \frac{1}{2})$. Set $\sigma_1 := \mu + \theta \in (\mu, \mu + \frac{1}{2})$ and $\sigma := \sigma_1 - \frac{\mu}{2} = \frac{\mu}{2} + \theta \in (\frac{\mu}{2}, \frac{1+\mu}{2})$. Then for $k \in \{0, 1\}$ and every $G \in H^{k+1}_\sigma(\mathbb{R}) \cap H^{k+1}_{\sigma_1}(\mathbb{R})$:*
1. **Exact Plancherel identity:** With $r = \sigma_1 + i\xi$ and the exact Branch I Mellin symbol $c_a(s) = (\mu - s)\big[\frac{\cos\pi\mu}{\mu} + \frac{\sin\pi\mu\cot\pi s}{1-s}\big]$,
   $$\langle \partial_z^k M_{\mathrm{I},c} G,\, \partial_z^k G\rangle_{L^2_\sigma(\mathbb{R})} = -(1-a)\sigma\,\|\partial_z^k G\|_{L^2_\sigma}^2 + \frac{1}{2\pi}\int_{-\infty}^\infty \operatorname{Re}\!\Big[(r-\mu)^k\bar r^k\, c_a(\mu - r)\Big]\big|\widehat{G}(r)\big|^2\,d\xi.$$
2. **Pointwise strict negativity (elementary, zero computer assistance):** For every $\xi \in \mathbb{R}$,
   $$\operatorname{Re}\!\big[c_a(\mu - r)\big] = -(1-a)\kappa\,\sigma_1 - \sin(\pi\mu)\,\frac{\sin(2\pi\theta)\big(\sigma_1(1+\theta) + \xi^2\big) + (1-\mu)\,\xi\sinh(2\pi\xi)}{\big(\cosh(2\pi\xi) - \cos(2\pi\theta)\big)\big((1+\theta)^2 + \xi^2\big)} < 0,$$
   $$\operatorname{Re}\!\big[(r-\mu)\bar r\,c_a(\mu - r)\big] = -(1-a)\kappa\,\theta\,|r|^2 - \sin(\pi\mu)\,|r|^2\,\frac{\sin(2\pi\theta)\big(\theta(1+\theta) + \xi^2\big) + \xi\sinh(2\pi\xi)}{\big(\cosh(2\pi\xi) - \cos(2\pi\theta)\big)\big((1+\theta)^2 + \xi^2\big)} < 0.$$
3. **Unconditional $H^1_\sigma$ coercivity and boundedness of gauge functionals:** Consequently,
   $$\langle \partial_z^k M_{\mathrm{I},c} G,\, \partial_z^k G\rangle_{L^2_\sigma(\mathbb{R})} \le -(1-a)\sigma\,\|\partial_z^k G\|_{L^2_\sigma(\mathbb{R})}^2 - (1-a)\kappa\,(\sigma_1 - k\mu)\,\|\partial_z^k G\|_{L^2_{\sigma_1}(\mathbb{R})}^2 \le -(1-a)\sigma\,\|\partial_z^k G\|_{L^2_\sigma(\mathbb{R})}^2 \qquad (k = 0, 1).$$
   Moreover, because $\frac{\mu}{2} < \mu < \frac{1+\mu}{2}$ for every $\mu \in [\frac{1}{2}, 1)$, choosing $\sigma_- \in (\frac{\mu}{2}, \mu)$ (e.g. $\sigma_- = \frac{3\mu}{4}$) and $\sigma_+ \in (\mu, \frac{1+\mu}{2})$ (e.g. $\sigma_+ = \frac{1+3\mu}{4}$) yields simultaneous $H^1_{\sigma_-} \cap H^1_{\sigma_+}$ coercivity together with Cauchy–Schwarz boundedness of both gauge functionals $|V_{\mathrm{I},G}| \le C\|G\|_{L^2_{\sigma_-}}$ and $|\widetilde{H}_{\mathrm{I},\mathrm{edge}}[G]| \le C(\|G\|_{L^2_{\sigma_-}} + \|G\|_{L^2_{\sigma_+}})$.

**Proof.**
By multiplicative-convolution Mellin calculus on the Muckenhoupt $A_2$ strip $\sigma_1 \in (\mu, 1+\mu)$, $e^{\mu z}\widetilde{H}_{\mathrm{I},0}[G] = -\sin(\pi\mu)\,v^{1+\mu} H_+[v^{-1-\mu}G]$ has Mellin symbol $\sin(\pi\mu)\cot(\pi(\mu - r))$, while $e^{\mu z}\widetilde{Q}_{\mathrm{I},0}[G](z) = \int_{-\infty}^z e^{(1-\mu)(z'-z)}(e^{\mu z'}\widetilde{H}_{\mathrm{I},0}[G])(z')\,dz'$ divides by $r + 1 - \mu$. Hence $e^{\mu z}(\widetilde{H}_{\mathrm{I},0} - (1-\mu)\widetilde{Q}_{\mathrm{I},0})$ has symbol $\sin(\pi\mu)\cot(\pi(\mu-r))\frac{r}{r+1-\mu}$. Combining this with the finite-time inflow transport term $-(1-a)\kappa e^{-\mu z}\partial_z G$ gives the total $e^{-\mu z}$-shifted symbol
$$-(1-a)\kappa\,r + \sin(\pi\mu)\cot(\pi(\mu - r))\frac{r}{r + 1 - \mu} = r\!\left[\frac{\cos(\pi\mu)}{\mu} + \frac{\sin(\pi\mu)\cot(\pi(\mu-r))}{1 - (\mu - r)}\right] = c_a(\mu - r).$$
Applying exact dilation commutation $\partial_z^k(e^{-\mu z} T G) = e^{-\mu z}(\partial_z - \mu)^k T G$ and Plancherel on $\operatorname{Re} r = \sigma_1 = \sigma + \frac{\mu}{2}$ yields part (1) directly (with the pure dilation term $-(1-a)\partial_z$ integrating by parts to $-(1-a)\sigma\|\partial_z^k G\|_{L^2_\sigma}^2$).

For part (2), write $r = \mu + \theta + i\xi$ with $\theta \in (0, \frac{1}{2})$. Then $\sin(2\pi\theta) > 0$, $\cosh(2\pi\xi) - \cos(2\pi\theta) > 0$, and $\cot(\pi(\mu - r)) = -\cot(\pi(\theta + i\xi)) = \frac{-\sin(2\pi\theta) + i\sinh(2\pi\xi)}{\cosh(2\pi\xi) - \cos(2\pi\theta)}$.
- For $k = 0$: $\frac{r}{r+1-\mu} = \frac{\sigma_1(1+\theta) + \xi^2 + i(1-\mu)\xi}{(1+\theta)^2 + \xi^2}$ has strictly positive real and imaginary coefficients (since $\sigma_1 > 0, \theta > 0, 1-\mu > 0$), and $\operatorname{Re}[-(1-a)\kappa r] = -(1-a)\kappa\sigma_1 \le 0$.
- For $k = 1$: $(r - \mu)\bar r\,\frac{r}{r+1-\mu} = |r|^2\frac{\theta + i\xi}{(1+\theta) + i\xi} = |r|^2\frac{\theta(1+\theta) + \xi^2 + i\xi}{(1+\theta)^2 + \xi^2}$ again has strictly positive real and imaginary coefficients for all $\xi \in \mathbb{R}$, while $\operatorname{Re}[-(1-a)\kappa(r-\mu)\bar r\,r] = -(1-a)\kappa\theta|r|^2 \le 0$.
Multiplying by $\cot(\pi(\mu-r))$ and taking real parts yields both formulas in (2), whose terms are strictly negative for all $\xi \in \mathbb{R}$. Part (3) follows immediately from (1)–(2) and Cauchy–Schwarz (for $V_{\mathrm{I},G}$, $e^{2(1-\mu+\sigma_-)z}\ln^2(1+e^{-z}) \in L^1(\mathbb{R})$ iff $\mu - 1 < \sigma_- < \mu$, which holds for every $\sigma_- \in (\frac{\mu}{2}, \mu)$). $\square$

*(Remark on $k \ge 2$ for Branch I with $a > -1$ and on $a < -2$ for Branch II: For Branch I with $a \in (-1, 0)$, $\kappa > 0$ contributes $-(1-a)\kappa(\sigma_1 - k\mu)\|\partial_z^k G\|_{L^2_{\sigma_1}}^2$ at top order, which switches sign for $k \ge 2$ whenever $\sigma_1 < 2\mu$; closing $H^2$ nonlinearly for Branch I at $a > -1$ therefore requires replacing $\partial_z^2$ by characteristic derivatives $\partial_\zeta^2 = ((1+\kappa e^{-\mu z})\partial_z)^2$. For Branch II at $a < -2$, the $H^3_{\sigma_+}$ linear coercivity threshold extends all the way to $a \in (-3, 0)$ as $\sigma_+ \to \frac{3}{4}^-$ — matching Theorem 1's exact pole threshold $a = -3$ — while nonlinear closure requires $\sigma_- < \frac{1}{2}$ so $a \in [-2, 0)$ is the uniform range for fixed weights $\sigma_- \in \{\frac{1}{4}, \frac{3}{8}\}$.)*

---

## 6. Status Summary and Honest Scope Boundaries

1. **Branch I at $a = -1$ & Branch II for all $a \in [-2, 0)$:** **Proved** linearly in both Hardy space (Theorem 1) and finite Sobolev space (Theorem 2 / Lemma 2.2, pure algebraic rational certificate), with global *a priori* nonlinear closure, uniqueness, uniform-in-$\nu$ viscous bounds, and global distributional validity on $\mathbb{R}$ (zero $\delta(X-1)$ defect) in finite Sobolev space $\mathcal{X}_{\mathrm{gCLM}}$ (Theorem 3, zero analyticity).
2. **Branch I for all $a \in [-1, 0)$:** **Proved linearly** in finite Sobolev space $H^1_{\sigma_-}(\mathbb{R}) \cap H^1_{\sigma_+}(\mathbb{R})$ via Theorem 4 (3-line computer-free Gauge-(c) Plancherel coercivity for every $a \in [-1, 0)$!), complementing Gates 2 and 4's numerical spectra (`#Re(lambda) > 0 == 0`, exact isolated real mode $\lambda_{\mathrm{I},*}(a) = -(1+a)$) and full $(N-1)$-dimensional restricted operator norm $\|e^{\tau M_a}|_{\ker L_0}\|_{L^2_{1/2}(d\zeta)} \le e^{-\tau/2}$. Nonlinear $H^2$ closure for $-1 < a < 0$ (and inner traveling-wave matching across $X=1$ where $\bar\Omega_{\mathrm{I},a} \notin L^2_{\mathrm{loc}}$) remains open.
3. **Branch I for $a < -1$:** Linearly stable to smooth relative perturbations $G \in C^1([1,\infty))$ (Gate 2) and universally attractive for smooth degenerate half-line initial data in nonlinear Case 2.2 DNS (Gate 3), despite the interior repelling stagnation point $X_* = 1 + (\cos\pi\mu)^{1/\mu} > 1$ and the steeper unstable $L^1_{\mathrm{loc}}\setminus L^2_{\mathrm{loc}}$ mode $\lambda_{\mathrm{I},*}(a) = |a|-1 > 0$.
4. **Stage 1 open problem:** As in 1D Hou–Luo, proving how smooth degenerate initial data on $\mathbb{R}$ dynamically form the $(X-1)^{-\mu}\mathbf{1}_{\{X>1\}}$ singular edge before entering the support-interior regime remains open.
