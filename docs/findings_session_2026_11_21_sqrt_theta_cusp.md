# Weak-Edge Compatibility via Transported $\sqrt{v}$ $\Theta$-Cusp (2026-11-21)

**Context:** Direct follow-up to Codex's 2026-11-21 scratch investigation (*"Weak-edge compatibility and modulation for the singular 1D Hou–Luo profile"*, Section 6), which proved:
- a gauge-independent obstruction $\frac{d}{d\tau}\log(B/A^2) = -Aa$ in the fixed half-power, pure-step/regular-$\Theta$ ansatz;
- edge-analytic linear rigidity ($G \equiv 0$) in the pure-step class;
- an infinite-dimensional non-analytic edge-flat compatible linear class ($a \equiv 0$) that fails nonlinearly via the strictly positive quadratic functional $J > 0$; and
- identified a square-root $\Theta$-cusp $\eta = D\sqrt{v} + o(\sqrt{v})$ (giving $(\log(B/A^2))' = -Aa - D/A$, requiring $D = -A^2 a$) as the sole surviving local mechanism (`[C]`), noting that a single fixed cusp shape fails.

**Evidence tiers:** `[D]` hand derivation under stated regularity; `[N]` 30–40 digit `mpmath` numerical check (`validate/edge_sqrt_theta_cusp_closure.py`, `TC-1..TC-3 PASS`). Not Lean-formalized; full nonlinear series convergence / Banach fixed-point existence remains open.

---

## 1. Exact Linear Closure for *Every* $G_0 \in \mathcal{X}$ (`[D, N]`)

Admit $\eta(v, \tau) = \sqrt{v}\,H(v, \tau)$ with $D(\tau) := H(0, \tau)$ and Gauge-(a) parameters $A \equiv B \equiv \xi \equiv 1$, $c_l = 2 - \varepsilon V_G$, $c_\omega = -c_l/2$ ($q = c_l + 2c_\omega \equiv 0$). At linear order:
$$\eta_\tau + 2v\eta_v = 0, \qquad G_\tau + 2vG_v = \tfrac{1}{2}\widetilde{Q}G + \bigl(\sqrt{v}\,\eta_v - \tfrac{1}{2}D\bigr), \qquad a(\tau) + D(\tau) = 0.$$
For **any** initial profile $G_0 \in \mathcal{X}$ (including edge-analytic data such as $G_0(v) = v/(b+v)$, with no flatness or zero-moment constraints), define
$$\eta_0(v) := -\frac{1}{2}\int_0^v y^{-1/2}(QG_0)(y)\,dy, \qquad G(v,\tau) := G_0(e^{-2\tau}v), \qquad \eta(v,\tau) := \eta_0(e^{-2\tau}v).$$
Then:
1. $\eta_\tau + 2v\eta_v \equiv 0$ identically.
2. Since $(QG_0)(0) = a_0$, $\eta_0(v) = -a_0\sqrt{v} + O(v^{3/2})$, so $D(\tau) = -e^{-\tau}a_0 = -a(\tau)$ for all $\tau \ge 0$ (exact weak-edge compatibility across $X = 1$).
3. By dilation covariance of $Q$, $\sqrt{v}\,\eta_v(v,\tau) = -\frac{1}{2}(QG)(v,\tau)$ cancels the nonlocal forcing pointwise, so $G_\tau + 2vG_v - \frac{1}{2}QG - \sqrt{v}\eta_v \equiv 0$.
4. The weighted norm decays strictly at the optimal convective rate: $\|G(\cdot,\tau)\|_{L^2_\sigma} = e^{-2\sigma\tau}\|G_0\|_{L^2_\sigma}$ for all $0 < \sigma < 1$.

**Check (`TC-1` on $G_0(v) = v/(2+v)$):** $\max|\Theta\text{ res}| = 2.8\times 10^{-17}$, $\max|\Omega\text{ res}| = 1.5\times 10^{-17}$, $\max|a(\tau)+D(\tau)| = 2.8\times 10^{-15}$.

---

## 2. Breaking the $O(\varepsilon^2)$ Quadratic Obstruction via Tricomi Inversion (`[D, N]`)

At $O(\varepsilon^2)$ around $(G_1, \eta_1) = (G_0(e^{-2\tau}v), \eta_0(e^{-2\tau}v))$ with $w := e^{-2\tau}v$:
- The quadratic source terms in the $\Omega$ and $\Theta$ equations are $e^{-\tau}S_G(w)$ and $e^{-\tau}S_\eta(w)$, where
  $$S_G(w) = -(QG_0 - V_{G_0})w G_0'(w) + \tfrac{1}{2}(QG_0 - a_0)G_0(w), \qquad S_\eta(w) = -(QG_0 - V_{G_0})w\eta_0'(w).$$
- The particular $\Theta$-solution $\eta_2^{\mathrm{part}}(v,\tau) = e^{-\tau}\widetilde{\eta}_2(w)$ has $\widetilde{\eta}_2(w) = -S_\eta(w) = -\frac{1}{2}\sqrt{w}\,QG_0(w)(QG_0(w) - V_{G_0})$, giving $D_2^{\mathrm{part}}(\tau) = D_2^{\mathrm{tilde}} e^{-2\tau}$ with $D_2^{\mathrm{tilde}} = -\frac{1}{2}a_0(a_0 - V_{G_0})$.
- Without a homogeneous $\Theta$-cusp correction, the $e^{-\tau}$ response $\widetilde{G}_2 = -S_G$ leaves a nonzero $O(\varepsilon^2)$ defect $(D_2 + a_2)(\tau) = +0.056852819440\,e^{-2\tau}$ for $G_0 = v/(2+v)$.
- **Exact $O(\varepsilon^2)$ cancellation:** Add a homogeneous transported $\Theta$-correction $\psi_2(e^{-2\tau}v)$ with $f_2(w) := \sqrt{w}\psi_2'(w)$ and expand $G_2(v,\tau) = \sum_{k\ge 1} e^{-k\tau} g_{2,k}(w)$. Setting $g_{2,1} = -S_G - f_2$, the $k=2$ equation $-2g_{2,2} = \frac{1}{2}Q g_{2,1} + \sqrt{w}\widetilde{\eta}_2'$ vanishes identically (terminating the series at $k=1$: $g_{2,k} \equiv 0$ for all $k \ge 2$) whenever
  $$Q(f_2) = -Q(S_G) + 2\sqrt{w}\,\widetilde{\eta}_2'.$$
  Writing $f_2 = -S_G + \tilde{f}_2$ reduces this to $Q(\tilde{f}_2) = R_2(w) := 2\sqrt{w}\,\widetilde{\eta}_2'(w)$, i.e. $H[w^{-1/2}\tilde{f}_2] = h_{\mathrm{rhs}}(w) := (w R_2(w))'$, which is inverted explicitly by the **one-sided Tricomi formula**
  $$\tilde{f}_2(w) = -\frac{w}{\pi}\operatorname{PV}\int_0^\infty y^{-1/2}\frac{h_{\mathrm{rhs}}(y) - h_{\mathrm{rhs}}(w)}{w - y}\,dy.$$
  Because $\tilde{f}_2(w) = O(w)$ at $w=0$, $\psi_2(w) = O(w^{3/2})$ has zero $\sqrt{w}$ cusp ($D(\psi_2) = 0$), while evaluating $Q(\tilde{f}_2)(0) = R_2(0) = D_2^{\mathrm{tilde}}$ forces $a(g_{2,1}) + D_2^{\mathrm{tilde}} = 0$ identically for all $\tau \ge 0$.

**Checks (`TC-2`, `TC-3`):** Tricomi inversion recovers $G_0$ to $2.2\times 10^{-26}$; $\tilde{f}_2(w)$ is smooth, real, $O(w)$ at $w=0$ ($\tilde{f}_2(0.01)/0.01 = +0.293549$), and satisfies $R_2(0) = D_2^{\mathrm{tilde}} = -0.0776932156$ (direct quadrature of $a(\tilde{f}_2)$ truncated at $w=400$ gives $-0.0776915005$, diff $1.7\times 10^{-6}$).

---

## 3. Scope & What Remains Open

- **Established (`[D, N]`):** Exact linear global weak-edge compatibility across $X = 1$ for every $G_0 \in \mathcal{X}$ via the $\sqrt{v}$ $\Theta$-cusp, and exact $O(\varepsilon^2)$ formal closure via Tricomi inversion (eliminating the pure-step quadratic obstruction $J > 0$).
- **Open:** Proving convergence of the formal $\varepsilon$-series or constructing the full nonlinear coupled $(G, \eta)$ compatibility manifold in a Banach space (and connecting to smooth degenerate initial data of Conjecture 2.4).
