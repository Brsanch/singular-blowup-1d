# Global Nonlinear Weak-Edge Compatible Manifold via Autonomous $s$-Clock, Exact Conjugation, and Strip-Analytic Contraction (2026-11-21)

**Goal:** Determine whether the global weak-edge compatibility condition $B(\tau) \equiv A(\tau)^2$ across $X = 1$ (equivalently $a(\tau) + D(\tau) \equiv 0$ in the $\sqrt{v}$ $\Theta$-cusp gauge) for the 1D Hou–Luo singular profile admits a **full nonlinear compatible manifold** via a strip-analytic Cauchy–Kovalevskaya / Nishida–Ovsyannikov fixed-point argument, or encounters a sharp obstruction.

**Verdict (updated 2026-11-22 after hostile review `findings_session_2026_11_21_weak_edge_hostile_review.md`):** **Autonomous $s$-clock reduction & numerical Picard convergence confirmed (`[D, N]`); full strip-analytic contraction NOT yet proved as stated (`[C]` / proof sketch).** Combining three exact structural reductions — (1) an autonomous finite-horizon $s$-clock that eliminates $V_G(\tau)$ identically, (2) the defect variable $W := (\frac{1}{2} + D)\widetilde H + \frac{1}{2} Q\widetilde G$, and (3) the exact conjugation $\frac{1}{2}Q(\mathcal{N}h) = -\mathcal{N}(u h)$ ($u = w^{-1/2}$) — turns the coupled nonlinear system into a forward-in-$s$ conjugated Duhamel system on $s \in [0, s_\infty)$ ($s_\infty = 1 + O(\varepsilon)$, covering all physical time $\tau \in [0, \infty)$). **Hostile-review corrections (A1, A2, A5, B1):**
- For the unmodified linear cusp ($h_0 \equiv 0$, $W(\cdot,0) \equiv 0$), the defect norm on the second Mellin line $\operatorname{Re} m = \sigma_2 > \frac{1}{2}$ blows up like $s^{1 - 2\sigma_2}$ as $s \to 0^+$, so the contraction in §3 is **false in the stated cone norm** without either (i) modifying the initial cusp by the $G_0$-dependent $O(\varepsilon^2)$ term $h_0(u) = \pi c_0(0)/(u+1)$ or (ii) using a time-weighted weakly singular Nishida norm.
- On $\operatorname{Re} m = \sigma_2 > \frac{1}{2}$, $S_g(t)$ is not uniformly bounded for $a \neq 0$ ($\|S_g(t)U\|_{\sigma_2} \sim |a(U)|t^{-(2\sigma_2 - 1)}$); closing the weakly singular Nishida scheme remains open.
- Compatibility $W(0,s) \equiv 0$ ($G(0,\tau) \equiv 0$) is an **inflow boundary condition** imposed by working in $L^2_{1/4}$ (not unique on $v > 0$, where $Q(1) = 0$ admits arbitrary edge-amplitude drifts).

**Evidence tiers:** `[D]` hand derivation for §1–§2 (algebraic reductions); `[C]` proof sketch for §3 (Nishida contraction has open gaps A1–A2); `[N]` numerical verification (`validate/edge_nonlinear_cusp_manifold.py`, `NM-1..NM-3 PASS`). Not Lean-formalized beyond the pointwise identities in `lean/SingularBlowup/WeakEdgeManifolds.lean`.

---

## 1. Exact Nonlinear System and Autonomous $s$-Clock Reduction (`[D, N]`)

In the $\sqrt{v}$ $\Theta$-cusp gauge $A \equiv B \equiv \xi \equiv 1$, $c_l(\tau) = 2 - V_G(\tau)$, $c_\omega(\tau) = -c_l(\tau)/2$ ($q = c_l + 2c_\omega \equiv 0$), substituting $\Omega(1+v,\tau) = v^{-1/2}(1+G(v,\tau))\mathbf{1}_{v>0}$ and $\Theta(1+v,\tau) = (\frac{\pi}{2}+\eta(v,\tau))\mathbf{1}_{v>0}$ directly into `(eq:dyn)` gives the exact nonlinear coupled PDE on $v > 0$:
$$\begin{cases}
G_\tau + \bigl(2 - V_G(\tau) + QG\bigr)\,v G_v - \dfrac{1}{2}(QG)(1 + G) = \sqrt{v}\,\eta_v, \\[6pt]
\eta_\tau + \bigl(2 - V_G(\tau) + QG\bigr)\,v \eta_v = 0,
\end{cases} \tag{$\star$}$$
with $G(0,\tau) = 0$ and $\eta(v,\tau) = D(\tau)\sqrt{v} + O(v^{3/2})$ as $v \to 0^+$ (for which the $v \to 0^+$ balance of the first equation is $a(\tau) + D(\tau) = 0$, matching $B(\tau) \equiv A(\tau)^2 = 1$).

> **Correction to earlier $O(\varepsilon^2)$ note (`[D, N]`, `NM-1`):** `findings_session_2026_11_21_sqrt_theta_cusp.md` (Section 2) wrote $+\frac{1}{2}(QG_0 - a_0)G_0$ (borrowing $-a_0$ from Gauge-(a), where $c_\omega = -1 + \frac{1}{2}(V_G - a)$) and assumed $V_{G_1}(\tau) = e^{-\tau} V_{G_0}$. In fact, in the $c_\omega = -c_l/2$ gauge the quadratic term is $+\frac{1}{2}(QG_0)G_0$, and for $G_0(w) = w/(b+w)$,
> $$V_{G_1}(\tau) = -\frac{e^\tau}{\pi}\int_0^\infty w^{-1/2}\ln\!\left(1 + \frac{e^{-2\tau}}{w}\right) G_0(w)\,dw = -2 + 2\sqrt{b}\,e^\tau\ln\!\left(1 + \frac{e^{-\tau}}{\sqrt{b}}\right) \neq e^{-\tau} V_{G_0}$$
> (verified to $5.84\times 10^{-22}$ in `NM-1`, with non-exponential gap $4.01\times 10^{-2}$ at $\tau = 0.5$). Both issues are eliminated simultaneously by the autonomous $s$-clock below.

### Autonomous finite-horizon $s$-clock (`[D]`)
Define the comoving dilation scale and finite-horizon clock:
$$\lambda(\tau) := \exp\!\left(-\int_0^\tau \bigl(2 - V_G(\tau')\bigr)\,d\tau'\right), \qquad s(\tau) := \int_0^\tau \sqrt{\lambda(\tau')}\,d\tau' \in [0, s_\infty),$$
and set $w := \lambda(\tau)v$, $\widetilde G(w, s) := G(v, \tau)$, and $\sqrt{w}\,\widetilde H(w, s) := \eta(v, \tau)$ (so $\widetilde D(s) = \widetilde H(0, s)$ and $\tilde a(s) = (Q\widetilde G)(0, s)$). By exact dilation covariance $Q[\widetilde G(\lambda\cdot)](v) = \sqrt{\lambda}\,(Q\widetilde G)(\lambda v)$, dividing $(\star)$ by $\sqrt{\lambda(\tau)} > 0$ **eliminates $V_G(\tau)$ identically**:
$$\begin{cases}
\widetilde G_s + (Q\widetilde G)\,D\widetilde G = \dfrac{1}{2}(Q\widetilde G)(1 + \widetilde G) + \left(\dfrac{1}{2} + D\right)\widetilde H, \\[8pt]
\widetilde H_s + (Q\widetilde G)\left(\dfrac{1}{2} + D\right)\widetilde H = 0,
\end{cases} \qquad D := w\partial_w = \partial_z, \quad s \in [0, s_\infty). \tag{$\dagger$}$$
At $w = 0$, the second equation gives the exact scalar edge law $\widetilde D'(s) + \frac{1}{2}\tilde a(s)\widetilde D(s) = 0$; under compatibility $\tilde a(s) + \widetilde D(s) \equiv 0$, this integrates in closed form to
$$\widetilde D(s) = -\tilde a(s) = \frac{-a_0}{1 + \frac{1}{2}a_0 s}.$$
Once $(\widetilde G, \widetilde H)$ is solved on $[0, S]$ with $S > 1$, the physical scale $r(s) := \sqrt{\lambda(s)}$ and physical time $\tau(s)$ are reconstructed **a posteriori** from
$$\frac{dr}{ds} = -1 + \frac{1}{2}V_{\widetilde G(\cdot, s)}(r^2), \qquad r(0) = 1, \qquad \tau(s) = \int_0^s \frac{ds'}{r(s')}.$$
For small solutions ($\|\widetilde G\| = O(\varepsilon)$), $|V_{\widetilde G}(r^2)| \le C\varepsilon r$, so $r(s)$ hits $0$ linearly at a unique finite horizon $s_\infty = 1 + O(\varepsilon) < S$, sending $\tau(s) \to +\infty$ as $s \to s_\infty^-$ with $\lambda(\tau) \sim c^2 e^{-2\tau}$ and optimal convective decay $\|G(\cdot,\tau)\|_{L^2_\sigma} = \lambda(\tau)^\sigma\|\widetilde G(\cdot, s(\tau))\|_{L^2_\sigma} \le C\varepsilon e^{-2\sigma\tau}$.

---

## 2. Why the Formal $O(\varepsilon^n)$ Tricomi Series Diverges — and How Conjugation Resums It (`[D, N]`)

Define the **weak-edge defect**
$$W(w, s) := \left(\frac{1}{2} + D\right)\widetilde H(w, s) + \frac{1}{2}(Q\widetilde G)(w, s).$$
Then $\widetilde G(0, s) \equiv 0 \iff W(0, s) \equiv 0$, and $(\dagger)$ is **exactly equivalent** to:
$$\begin{cases}
\widetilde G_s = W + \mathcal{P}(\widetilde G), & \mathcal{P}(\widetilde G) := -(Q\widetilde G)\,D\widetilde G + \dfrac{1}{2}(Q\widetilde G)\,\widetilde G, \\[6pt]
W_s - \dfrac{1}{2} Q W = -\left(\dfrac{1}{2} + D\right)\bigl((Q\widetilde G) W\bigr) + \mathcal{R}(\widetilde G), & \mathcal{R}(\widetilde G) := \dfrac{1}{2}\left(\dfrac{1}{2}+D\right)\bigl((Q\widetilde G)^2\bigr) + \dfrac{1}{2}Q(\mathcal{P}(\widetilde G)).
\end{cases} \tag{$\ddagger$}$$
(Proof of equivalence: defining $\widetilde H$ by $\widetilde H_s = -(Q\widetilde G)(W - \frac{1}{2}Q\widetilde G)$ with $\widetilde H(\cdot,0) = (\frac{1}{2}+D)^{-1}(-\frac{1}{2}Q\widetilde G_0 + W_0)$, the residual $E := (\frac{1}{2}+D)\widetilde H + \frac{1}{2}Q\widetilde G - W$ satisfies $E(\cdot,0) = 0$ and $E_s \equiv 0$.)

### $O(\varepsilon^3)$ formal polynomial hierarchy & divergence diagnosis (`[D, N]`, `NM-1`)
Expanding $\widetilde G = \sum_{n\ge 1}\varepsilon^n \hat G_{n-1}$, $W = \sum_{n\ge 2}\varepsilon^n \hat W_{n-1}$ in polynomials in $s$ solves $W_{n,s} - \frac{1}{2}Q W_n = F_n(w,s)$ backward from the top power of $s$ via the formal Neumann series
$$W_n(w, s) = -\sum_{m=0}^{n-2} 2^{m+1}\,(Q^{-1})^{m+1}\bigl(\partial_s^m F_n(w, s)\bigr).$$
At $O(\varepsilon^3)$ ($n=3$), the polynomial hierarchy terminates cleanly at degree $2$ in $s$ (`NM-1`: $R_{2,2}(0) = -\frac{1}{4}a_0^3$ exact to $10^{-30}$, $g_{2,2} = -Q^{-1}(R_{2,2}) = O(w)$). However, stepping down from $k = n$ to $k = 1$ applies $Q^{-1}$ (a 1st-order operator in $z = \ln w$) $O(2n)$ times to $G_0$ while $s$-integration divides by only $n!$, producing Gevrey-1 factorial growth $\sim (2n)!/n! \sim 4^n n!$.

### Exact conjugation resums the divergent Tricomi series (`[D, N]`, `NM-2`)
Under the exact Mellin conjugation $W = \mathcal{N} h$ with multiplier $N(m) = \frac{\cos(\pi m)\Gamma(2m+2)}{\pi}$ ($u := w^{-1/2} = e^{-z/2}$), we have $\frac{1}{2}Q(\mathcal{N}h) = -\mathcal{N}(u h)$, so $W_s - \frac{1}{2}Q W = F$ becomes the **pointwise ODE**
$$h_s(u, s) + u\,h(u, s) = f(u, s), \qquad f := \mathcal{N}^{-1} F,$$
with exact forward Duhamel solution
$$h(u, s) = e^{-u s} h_0(u) + \int_0^s e^{-u(s - s')} f(u, s')\,ds', \qquad W(s) = S(s) W_0 + \int_0^s S(s - s') F(s')\,ds', \qquad S(s) := \mathcal{N} e^{-u s} \mathcal{N}^{-1}.$$
The divergent backward Tricomi series $-\sum 2^{m+1}(Q^{-1})^{m+1}\partial_s^m F$ is simply the **asymptotic $u \to \infty$ integration-by-parts expansion** of the convergent integral $\int_0^s e^{-u(s-s')} f(u, s')\,ds' \sim \sum_{m\ge 0}(-1)^m u^{-(m+1)}\partial_s^m f(u,s)$!

Crucially, even with **zero initial defect $h_0 \equiv 0$ ($W_0 \equiv 0$)**:
- For $F$ in R-class ($F(0, s) \neq 0$), $f(u, s) \to \pi F(0, s)$ as $u \to \infty$, so $h(u, s) \sim \frac{\pi F(0, s)}{u}$ as $u \to \infty$ for all $s > 0$.
- The simple pole of $\widehat h(m, s)$ at $m = \frac{1}{2}$ is cancelled by the simple zero $N(\frac{1}{2}) = 0$, leaving $W(w, s) = O(w)$ ($W(0, s) \equiv 0$) and $a(W)(s) = -\frac{2}{\pi}\lim_{u\to\infty} u\,h(u, s) = -2 F(0, s)$ **automatically** for all $s > 0$! (The $s = 0$ jump in $u h(u, s)$ is removable if desired by setting $h_0(u) = f(u, 0)/(u+1)$, which adds only a smooth $O(\varepsilon^2 w^{3/2})$ initial cusp with zero $\sqrt{w}$ component.)

---

## 3. Strip-Analytic Cauchy–Kovalevskaya / Nishida–Ovsyannikov Fixed-Point Proof Sketch (`[D]`)

Fix $\rho_0 \in (0, \pi)$ (e.g. $\rho_0 = 1$) and define the strip-analytic intersection space $\mathcal{X}_\rho := H^2_\rho(\frac{1}{4}) \cap H^2_\rho(\frac{3}{4})$ on $\mathbb{S}_\rho := \{|\operatorname{Im} z| < \rho\}$, with norm $\|f\|_\rho^2 := \int_{\mathbb{R}} e^{2\rho|\xi|}(|\widehat f(\frac{1}{4}+i\xi)|^2 + |\widehat f(\frac{3}{4}+i\xi)|^2)\,d\xi$.

1. **Algebra & $Q$-smoothing (`[D]`):** Since $\mathcal{M}[Q G](\zeta) = \frac{\cot(\pi\zeta)}{\zeta+1}\widehat G(\zeta+\frac{1}{2})$ gains 1 derivative in $\xi$, decomposing $QG = a(G) + \widetilde Q G$ (with $|a(G)| \le C\|G\|_\rho$ and $\widetilde Q G \in H^2_\rho(\frac{1}{4})$ after crossing the simple pole of $\cot(\pi\zeta)$ at $\zeta = 0$, while $QG \in H^2_\rho(-\frac{1}{4})$) yields $\|Q G\|_{H^\infty(\mathbb{S}_\rho)} \le C\|G\|_\rho$ with **zero radius shrinkage**, so multiplication by $Q G$ is bounded on $\mathcal{X}_\rho$ and $\mathcal{R}_\rho := H^2_\rho(-\frac{1}{4})$.
2. **Strip-analytic bounds on $S(s)$ (`[D, N]`, `NM-2`; open on $\sigma_2 > \frac{1}{2}$ per A1–A2):** On $\mathbb{S}_\rho$ ($\rho < \pi$), $\operatorname{Re}(u) = e^{-\operatorname{Re}(z)/2}\cos(\operatorname{Im}(z)/2) \ge \cos(\rho/2)|u| > 0$. Hence for all $s \ge 0$ and all $k \ge 0$:
   $$|e^{-u s}| \le 1, \qquad |\partial_z^k(e^{-u s})| \le C_k(\rho), \qquad \int_0^s |\partial_z^k(u e^{-u(s-s')})|\,ds' \le \widetilde C_k(\rho).$$
   By Stirling's formula, $|N(\sigma + i\xi)| \sim C|\xi|^{2\sigma + 3/2}$ (the $e^{\pm\pi|\xi|}$ factors of $\cos$ and $\Gamma$ cancel identically). On the first line $\sigma_1 = \frac{1}{4}$, $S_g(s)$ has 0 derivative loss uniformly in $s \ge 0$ and $\mathcal{T}$ loses 1 derivative. **Caveat (A1–A2):** On the second line $\sigma_2 > \frac{1}{2}$, $\|S_g(t)U\|_{\sigma_2} \sim |a(U)|t^{-(2\sigma_2 - 1)}$ and $\mathcal{T}$ blows up like $s^{1-2\sigma_2}$ as $s \to 0^+$ when $h_0 \equiv 0$; closing the contraction on $\mathcal{X}_\rho$ requires either the $O(\varepsilon^2)$ initial cusp modification $h_0(u) = \pi c_0(0)/(u+1)$ or a time-weighted weakly singular Nishida norm (not yet written).
3. **Weighted cone contraction (`[D]`):** Writing $\widetilde G = \varepsilon\hat G$, $W = \varepsilon^2\hat W$, the Duhamel map $(\clubsuit)$ is
   $$\hat G(s) = G_0 + \varepsilon\int_0^s \bigl[\hat W(s') + \mathcal{P}(\hat G(s'))\bigr]\,ds', \qquad \hat W(s) = -\varepsilon\int_0^s S_g(s-s')\left(\frac{1}{2}+D\right)\bigl((Q\hat G)\hat W\bigr)\,ds' + \mathcal{T}\bigl[\mathcal{R}(\hat G)\bigr](s).$$
   On the cone $0 \le s < \rho_0/a$, $0 < \rho < \rho_0 - a s$ with $d(\rho, s) := \rho_0 - \rho - a s$, endow $(\delta\hat G, \delta\hat W)$ with the Nishida–Ovsyannikov norm
   $$\|(\delta\hat G, \delta\hat W)\|_* := \sup_{s,\rho}\Bigl[d(\rho,s)^{1/2}\|\delta\hat G(s)\|_\rho + \gamma\,d(\rho,s)^{3/2}\|\delta\hat W(s)\|_\rho\Bigr].$$
   Assuming $\|\hat G(s)\|_\rho \le M_G$ and $d(\rho,s)\|\hat W(s)\|_\rho \le M_W$ (maintained across Nishida's $a_j \uparrow a_\infty$ sequence with $a_{j+1} - a_j = c/(j+1)^2$), choosing $\rho'(s') = \rho + \frac{1}{2}d(\rho, s')$ yields the contraction factor
   $$\kappa \le \frac{2\varepsilon}{\gamma a} + 2\sqrt{2}\,\gamma\,C_T M_G + O\!\left(\frac{\varepsilon}{a}\right).$$
   Setting $\gamma := (20 C_T M_G)^{-1}$ and $a := 40\varepsilon/\gamma = 800 C_T M_G \varepsilon$ gives $\kappa < \frac{1}{2}$ on $[0, \rho_0/a)$. For $\varepsilon < \varepsilon_0(\rho_0, \|G_0\|_{\rho_0})$, the existence horizon $T = \rho_0/a = \frac{\rho_0}{800 C_T M_G \varepsilon} > 2 > s_\infty = 1 + O(\varepsilon)$, yielding a unique strip-analytic global solution on all of $\tau \in [0, \infty)$.

---

## 4. Verification Summary (`validate/edge_nonlinear_cusp_manifold.py`, `NM-1..NM-3`)

| Check | Content | Result |
|---|---|---|
| `NM-1` | 35-digit `mpmath`: exact $V_{G_1}(\tau)$ vs quadrature (non-exponential gap $4.01\times 10^{-2}$); $O(\varepsilon^2, \varepsilon^3)$ edge identities & $g_{2,2} = O(w)$ via Tricomi | `PASS`: $V_{G_1}$ err `5.84e-22`, edge err `0.0e+00`, $g_{2,2}(0.001)/0.001 = +0.402676$ |
| `NM-2` | Conjugation $\frac{1}{2}Q(\mathcal{N}h) = -\mathcal{N}(uh)$; strip-analytic kernel bounds on $|\operatorname{Im}z| \le 1$; zero-initial-data Duhamel on $N_2$ | `PASS`: `conj_err = 2.69e-9`, kernel bounds sharp, `N2_PDE_res = 7.61e-6`, `|0.5*QW2(0)+N2(0)| = 3.03e-7` |
| `NM-3` | Full nonlinear conjugated Picard iteration ($\varepsilon = 0.05$, $G_0 = w/(2+w)^2$, $W(\cdot,0) \equiv 0$, $s \in [0, 1.2] > s_\infty = 0.998662$) + physical-time $(\star)$ check across $\tau \in \{0.2, 0.6, 1.5, 3.0\}$ | `PASS`: `||dG|| = 1.24e-12`, `||dW|| = 3.98e-12`; $s$-clock `res_G = 9.93e-6`, `res_H = 3.89e-9`, `|a+D| = 1.72e-7`, `D_end_err = 1.91e-9`; physical-time `res_G(*) = 2.61e-6`, `res_eta(*) = 9.16e-6`, `|a(tau)+D(tau)| = 1.35e-7` |

---

## 5. Scope & What Remains Open

- **Established (`[D, N]`):** Exact autonomous finite-horizon $s$-clock reduction $(\star) \iff (\dagger)$, exact defect/conjugation resummation of the divergent formal Tricomi series, and numerical Picard convergence on clean domains `[-64, 64]` with $V_G$-sensitive physical-time verification.
- **Open (per hostile review A1–A2, A5, B3):** (1) Closing the strip-analytic contraction on the second line $\sigma_2 > \frac{1}{2}$ (either via modified initial cusp $h_0(u) = \pi c_0(0)/(u+1)$ or a time-weighted weakly singular Nishida norm); (2) uniqueness of weak solutions under the inflow boundary condition $G(0,\tau) = 0$; (3) verifying the distributional singular-product convention against arXiv:2604.01868 and connecting to smooth degenerate initial data of Conjecture 2.4.
