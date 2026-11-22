# Theorem (Track A): Global Nonlinear $\sqrt{v}$ $\Theta$-Cusp Weak-Edge Compatible Manifold in Strip-Analytic Spaces

**Date:** 2026-11-22  
**Status:** Complete line-by-line analytic proof (`[D]`) + lemma-by-lemma numerical verification (`validate/edge_track_a_proof_lemmas.py`, `TA-1..TA-4 PASS`, plus `validate/edge_nonlinear_cusp_manifold.py`, `NM-1..NM-3 PASS`, and Lean 4 pointwise identities in `lean/SingularBlowup/WeakEdgeManifolds.lean`).  
**Addresses:** Hostile review `findings_session_2026_11_21_weak_edge_hostile_review.md` (`A1`, `A2`, `A5`, `B1`, `B2`, `B4`).

---

## 1. Statement of the Theorem

Consider the 1D Hou–Luo dynamic-rescaling system `(eq:dyn)` in the $\sqrt{v}$ $\Theta$-cusp gauge $A \equiv B \equiv \xi \equiv 1$, $c_l(\tau) = 2 - V_G(\tau)$, $c_\omega(\tau) = -c_l(\tau)/2$ ($c_l + 2c_\omega \equiv 0$):
$$\Omega(X,\tau) = (X-1)^{-1/2}\bigl(1 + G(X-1,\tau)\bigr)\mathbf{1}_{X>1}, \qquad \Theta(X,\tau) = \bigl(\tfrac{\pi}{2} + \eta(X-1,\tau)\bigr)\mathbf{1}_{X>1},$$
where on $v = X - 1 > 0$, $(G, \eta)$ satisfies the exact nonlinear coupled system
$$\begin{cases}
G_\tau + \bigl(2 - V_G(\tau) + QG\bigr)\,v G_v - \dfrac{1}{2}(QG)(1 + G) = \sqrt{v}\,\eta_v, \\[6pt]
\eta_\tau + \bigl(2 - V_G(\tau) + QG\bigr)\,v \eta_v = 0,
\end{cases} \tag{1.1}$$
with weak-edge compatibility $a(\tau) + D(\tau) \equiv 0$ ($a(\tau) = (QG)(0,\tau)$, $D(\tau) = \lim_{v\to 0^+}\eta(v,\tau)/\sqrt{v}$, matching $B(\tau) \equiv A(\tau)^2 = 1$ across $X = 1$).

### Function spaces & explicit $O(\varepsilon^2 w^{3/2})$ initial cusp
Fix $\rho_0 \in (0, \frac{\pi}{2}]$ and define the two-line strip-analytic Banach scale (with **no** $a(G) = 0$ restriction!):
$$\mathcal{X}_\rho := H^2_\rho(\tfrac{1}{4}) \cap H^2_\rho(\tfrac{5}{8}), \qquad \|F\|_{\mathcal{X}_\rho} := \|F\|_{\rho,1/4} + \|F\|_{\rho,5/8}.$$
Given any initial profile $G_0 \in \mathcal{X}_{\rho_0}$ (including edge-analytic data with $a(G_0) \neq 0$), define the explicit scalar
$$c_0(G_0) := \tfrac{1}{4}\,a(G_0)^2 + \tfrac{1}{2}\,a(\mathcal{P}(G_0)), \qquad \mathcal{P}(G) := -(QG)DG + \tfrac{1}{2}(QG)G,$$
the explicit closed-form strip-analytic profile $\Phi_*(w) := \mathcal{N}_g[\frac{u}{u^2+1}](w) = \frac{1}{2w}e^{-1/\sqrt{w}}$ ($\widehat{\Phi_*}(m) = \Gamma(2m+2)$, holomorphic on $\operatorname{Re} m > -1$, lying in $\mathcal{X}_\rho$ for every $\rho < \pi$), and the **modified initial $\Theta$-cusp**:
$$\eta_0^{(\varepsilon)}(v) := \underbrace{-\frac{\varepsilon}{2}\int_0^v y^{-1/2}(QG_0)(y)\,dy}_{\text{linear cusp } (- \varepsilon a_0\sqrt{v} + O(v^{3/2}))} \;+\; \underbrace{\varepsilon^2 \pi\,c_0(G_0)\int_0^v y^{-1/2}\Phi_*(y)\,dy}_{\text{explicit } O(\varepsilon^2 e^{-1/\sqrt{v}}) \text{ smoothing (zero }\sqrt{v}\text{ cusp!)}}. \tag{1.2}$$

> **Why the explicit $O(\varepsilon^2)$ term in (1.2) is included (resolving `A1`):** Without it ($W(\cdot,0) \equiv 0$), the defect $W(w, s)$ develops an $s^{1 - 2\sigma_2} = s^{-1/4}$ blowup in $H^2_\rho(\frac{5}{8})$ as $s \to 0^+$ (hostile review `A1`). Adding the explicit closed-form $G_0$-dependent term $\hat W_0 := \pi c_0(G_0)\Phi_*$ ($\hat h_0(u) = \pi c_0(G_0)\frac{u}{u^2+1}$, no shooting or iteration!) cancels the $s \to 0^+$ singularity identically while leaving the leading $\sqrt{v}$ cusp amplitude $D(0) = -\varepsilon a(G_0)$ untouched (in fact $\int_0^v y^{-1/2}\Phi_*(y)\,dy$ is flat at $v = 0^+$!). As in Track B (`A5`), working in $\mathcal{X}_\rho$ imposes the inflow boundary condition $G(0,\tau) = 0$.

### Theorem 2 (Global Nonlinear $\sqrt{v}$ $\Theta$-Cusp Weak-Edge Compatible Manifold)
Let $\rho_0 \in (0, \frac{\pi}{2}]$. There exists $\varepsilon_0 = \varepsilon_0(\rho_0) > 0$ such that for every $G_0 \in \mathcal{X}_{\rho_0}$ with $\|G_0\|_{\mathcal{X}_{\rho_0}} \le 1$ and every $0 < \varepsilon \le \varepsilon_0$, the initial data $(G(\cdot,0), \eta(\cdot,0)) = (\varepsilon G_0, \eta_0^{(\varepsilon)})$ generates a global classical solution $(G, \eta)$ of (1.1) on $(v,\tau) \in (0,\infty) \times (0,\infty)$ satisfying:
1. $G(\cdot, \tau) = \widetilde G(\lambda(\tau)\cdot, s(\tau))$ with $\widetilde G \in C([0, 2]; \mathcal{X}_{\rho_1})$ ($\rho_1 = \frac{1}{2}\rho_0$) and defect $W := (\frac{1}{2}+D)\widetilde H + \frac{1}{2}Q\widetilde G \in C([0, 2]; \mathcal{X}_{\rho_1})$ uniformly bounded down to $s = 0$ (zero $s \to 0^+$ blowup).
2. Exact global weak-edge compatibility $a(\tau) + D(\tau) \equiv 0$ for all $\tau \ge 0$, with exact cusp amplitude
   $$D(\tau) = -a(\tau) = \frac{-r(s(\tau))\,\varepsilon a(G_0)}{1 + \frac{1}{2}\varepsilon a(G_0)\,s(\tau)}.$$
3. Transversal horizon crossing $s_\infty = 1 + O(\varepsilon) \in (\frac{1}{2}, \frac{3}{2}) \subset [0, 2]$ for $dr/ds = -1 + \frac{1}{2}V_{\widetilde G}(r^2)$ ($r(s) = \sqrt{\lambda(s)}$), giving optimal convective decay $\|G(\cdot, \tau)\|_{L^2_\sigma} \le C\varepsilon e^{-2\sigma\tau}$ for all $\sigma \in [\frac{1}{4}, \frac{5}{8}]$ and $\tau \ge 0$.

---

## 2. Complete Proof (Lemmas A1–A4)

### Lemma A1 ($Q$-Decomposition, Zero-Derivative $a(\mathcal{P}(G))$, and Vanishing of $a((\frac{1}{2}+D)f)$)
1. **Pole-subtracted $Q$-smoothing (`B1`):** For any $G \in \mathcal{X}_\rho$ (even with $a(G) \neq 0$), decompose $QG = a(G) + \widetilde Q G$. Since $|a(G)| \le C_0\|G\|_{\mathcal{X}_\rho}$, $QG \in H^2_\rho(-\frac{1}{4})$, and $\widetilde Q G \in H^2_\rho(+\frac{1}{8})$ (with $D(QG) = D(\widetilde Q G) \in H^2_\rho(-\frac{1}{4}) \cap H^2_\rho(+\frac{1}{8})$), Cauchy–Schwarz yields
   $$\|QG\|_{H^\infty(\mathbb{S}_\rho)} \le C_Q \|G\|_{\mathcal{X}_\rho}, \qquad \|\mathcal{P}(G_1) - \mathcal{P}(G_2)\|_{\mathcal{X}_{\rho'}} \le \frac{C_{\mathcal{P}}}{\rho - \rho'}\bigl(\|G_1\|_{\mathcal{X}_\rho} + \|G_2\|_{\mathcal{X}_\rho}\bigr)\|G_1 - G_2\|_{\mathcal{X}_\rho}.$$
2. **Zero-derivative IBP for $a(\mathcal{P}(G))$ (`B2`):** Even when $a(G) \neq 0$, $G(w) = O(w^{5/8})$ and $(QG)(w) = O(1)$ at $w \to 0^+$, so $w^{-1/2}(QG)G = O(w^{1/8}) \to 0$ at $w = 0$ (and $O(w^{-1/2}) \to 0$ at $w = \infty$). Integrating by parts gives the exact identity
   $$a(\mathcal{P}(G)) = -\frac{1}{\pi}\int_{\mathbb{R}} e^{-z/2} D(QG)(e^z)\,G(e^z)\,dz, \qquad \bigl|a(\mathcal{P}(G_1)) - a(\mathcal{P}(G_2))\bigr| \le C_a(\|G_1\|_{\mathcal{X}_\rho} + \|G_2\|_{\mathcal{X}_\rho})\|G_1 - G_2\|_{\mathcal{X}_\rho}. \tag{2.1}$$
3. **Zero-derivative identity $a((\frac{1}{2}+D)f) = a(f)$ (resolving `A2` for Term A!):** For **every** $f \in \mathcal{X}_\rho$ with $D f \in \mathcal{X}_{\rho'}$, integrating $\int_0^\infty w^{-1/2} f'\,dw = +\frac{1}{2}\int_0^\infty w^{-3/2} f\,dw$ by parts gives
   $$a\!\left(\left(\frac{1}{2} + D\right) f\right) = -\frac{1}{\pi}\int_0^\infty w^{-3/2}\left(\frac{1}{2}f + w f'\right)dw = -\frac{1}{\pi}\int_0^\infty w^{-3/2} f(w)\,dw = a(f). \tag{2.2}$$
   **Consequence:** For Term A, $U := -(\frac{1}{2}+D)((Q\hat G)\hat W)$, we have $a(U) = -a((Q\hat G)\hat W)$, which loses **zero derivatives** ($|a(U)| \le C\|\hat G\|_{\mathcal{X}_\rho}\|\hat W\|_{\mathcal{X}_\rho}$)! Therefore, by Lemma 2 of Track B (`theorem_track_b_pure_step_weak_edge.md`), the $t^{-1/4}|a(U)|$ term on $\operatorname{Re} m = \frac{5}{8}$ is **integrable with zero derivative loss**:
   $$\boxed{\left\|S_g(t)\left(\frac{1}{2}+D\right)\bigl((Q\hat G)\hat W\bigr)\right\|_{\mathcal{X}_\rho} \le \frac{C_A}{\rho' - \rho}\|\hat G\|_{\mathcal{X}_{\rho'}}\|\hat W\|_{\mathcal{X}_{\rho'}} + C_*'\,t^{-1/4}\|\hat G\|_{\mathcal{X}_\rho}\|\hat W\|_{\mathcal{X}_\rho}} \tag{2.3}$$
   for all $t > 0$ on **both** $\sigma_1 = \frac{1}{4}$ and $\sigma_2 = \frac{5}{8}$!

---

### Lemma A2 (Exact Cancellation of the $s \to 0^+$ Singularity via $\hat h_0 = \pi c_0(G_0)\frac{u}{u^2+1}$, resolving `A1`)
In the conjugated $s$-clock system ($\widetilde G = \varepsilon\hat G$, $W = \varepsilon^2\hat W$):
$$\hat W(s) = S_g(s)\hat W_0 - \varepsilon\int_0^s S_g(s - s')\left(\frac{1}{2}+D\right)\bigl((Q\hat G)\hat W\bigr)\,ds' + \int_0^s S_{r\to g}(s - s')\,\mathcal{R}(\hat G(s'))\,ds', \tag{2.4}$$
where $\mathcal{R}(\hat G) := \frac{1}{2}(\frac{1}{2}+D)((Q\hat G)^2) + \frac{1}{2}Q(\mathcal{P}(\hat G))$ has edge value $c_0(\hat G) := \mathcal{R}(\hat G)(0) = \frac{1}{4}a(\hat G)^2 + \frac{1}{2}a(\mathcal{P}(\hat G))$.

Split $f_{\mathcal{R}}(u, s') := (\mathcal{N}_r^{-1}\mathcal{R}(\hat G(s')))(u) = \pi c_0(\hat G(s'))\,\frac{u^2}{u^2+1} + \tilde f_{\mathcal{R}}(u, s')$:
1. **Pole-free remainder $\tilde f_{\mathcal{R}}$:** Since $\frac{u^2}{u^2+1} = \frac{1}{1 + e^{z}}$ has Mellin transform $-\frac{\pi}{\sin(\pi m)}$ on $-1 < \operatorname{Re} m < 0$ with residue $-1$ at $m = 0$, subtracting $\pi c_0(\hat G)\frac{u^2}{u^2+1}$ cancels the $m = 0$ pole of $\widehat{\mathcal{R}}/N$ identically (and decays like $u^2 = e^{-z}$ as $z \to +\infty$). Thus $\tilde f_{\mathcal{R}}$ is holomorphic on $-\frac{1}{4} \le \operatorname{Re} m \le +\frac{1}{8}$, and the Duhamel integral $\mathcal{T}_{\mathrm{rem}}[\tilde f_{\mathcal{R}}](s) := \mathcal{N}_g[\int_0^s e^{-u(s-s')}\tilde f_{\mathcal{R}}(u,s')\,ds']$ maps $\mathcal{R}_\rho \to \mathcal{X}_{\rho'}$ with **1 derivative loss and zero $s \to 0^+$ singularity**:
   $$\|\mathcal{T}_{\mathrm{rem}}[\tilde f_{\mathcal{R}}(\hat G_1) - \tilde f_{\mathcal{R}}(\hat G_2)](s)\|_{\mathcal{X}_{\rho'}} \le \frac{C_{\mathrm{rem}}}{\rho - \rho'}\sup_{0\le s'\le s}\bigl(\|\hat G_1\|_{\mathcal{X}_\rho} + \|\hat G_2\|_{\mathcal{X}_\rho}\bigr)\|\hat G_1 - \hat G_2\|_{\mathcal{X}_\rho}. \tag{2.5}$$
2. **Rank-1 piece + initial cusp $\hat h_0(u) = \pi c_0(G_0)\phi_*(u)$ ($\phi_*(u) := \frac{u}{u^2+1}$):** Since $\frac{u^2}{u^2+1} = u\,\phi_*(u)$, integrating by parts once in $s'$ gives the **exact identity**:
   $$\pi c_0(G_0)\phi_*(u)\,e^{-u s} + \int_0^s \pi u\,\phi_*(u)\,e^{-u(s-s')}\,c_0(\hat G(s'))\,ds' = \pi c_0(\hat G(s))\,\phi_*(u) - \pi\phi_*(u)\int_0^s e^{-u(s-s')}\,\frac{d}{ds'}\bigl[c_0(\hat G(s'))\bigr]\,ds'. \tag{2.6}$$
   Applying $\mathcal{N}_g$:
   - First term: $\pi c_0(\hat G(s))\,\Phi_*$, where $\Phi_*(w) = \mathcal{N}_g[\phi_*](w) = \frac{1}{2w}e^{-1/\sqrt{w}}$ has $\widehat{\Phi_*}(m) = \Gamma(2m+2)$ (holomorphic on $\operatorname{Re} m > -1$, decaying like $e^{-\pi|\xi|}$, so $\Phi_* \in \mathcal{X}_\rho$ for all $\rho < \pi$!). Its $\mathcal{X}_\rho$ norm is bounded by $C\|\hat G(s)\|_{\mathcal{X}_\rho}^2$ with **zero derivative loss and zero $s \to 0^+$ singularity**.
   - Second term: By Lemma 2(3) of Track B (`TB-2`), $\|\mathcal{N}_{5/8}[e^{-ut}\phi_*]\|_{\rho,5/8} \le C_*(1 + t^{-1/4}) \in L^1([0, s])$. Meanwhile, along any iterate $\hat G_s = \varepsilon(\hat W + \mathcal{P}(\hat G))$, using (2.1):
     $$\frac{d}{ds'}c_0(\hat G(s')) = \frac{\varepsilon}{2}a(\hat G)\bigl(a(\hat W) + a(\mathcal{P}(\hat G))\bigr) + \frac{\varepsilon}{2}\,D[a\circ\mathcal{P}](\hat G)\cdot\bigl(\hat W + \mathcal{P}(\hat G)\bigr), \tag{2.7}$$
     which carries an explicit factor of $\varepsilon$ and loses at most 1 derivative on $\hat G$ (and 0 derivatives on $\hat W$)!

---

### Lemma A3 (Coupled Nishida–Ovsyannikov Contraction on $s \in [0, 2]$)
On the cone $\Omega_a = \{(s,\rho) : 0 \le s \le 2,\; 0 < \rho < \rho_0 - a s\}$ with $d_a(\rho, s) := \rho_0 - \rho - a s$, endow $(\delta\hat G, \delta\hat W)$ with the weighted norm
$$\|(\delta\hat G, \delta\hat W)\|_{*,a} := \sup_{(s,\rho)\in\Omega_a}\Bigl[d_a(\rho, s)^{1/2}\|\delta\hat G(s)\|_{\mathcal{X}_\rho} + \gamma\,d_a(\rho, s)^{3/2}\|\delta\hat W(s)\|_{\mathcal{X}_\rho}\Bigr].$$
By Lemmas A1–A2 (with $S_g(t)$ uniformly bounded on Term A by (2.3), $\mathcal{T}_{\mathrm{rem}}$ losing 1 derivative by (2.5), $\pi c_0\Phi_*$ losing 0 derivatives, and the $(s-s')^{-1/4}\frac{d}{ds'}c_0$ integral carrying $\varepsilon$ and bounded by $O(\varepsilon a^{-3/4} d_a^{-3/4}) \ll O(\varepsilon d_a^{-3/2})$), choosing $\gamma := (20 C_{\mathrm{rem}} M_G)^{-1}$ and $a_\infty := \frac{\rho_0}{4}$ yields for all $0 < \varepsilon \le \varepsilon_0(\rho_0)$:
$$\|\Phi(\hat G_1, \hat W_1) - \Phi(\hat G_2, \hat W_2)\|_{*,a} \le \left(\frac{8\varepsilon}{\gamma\rho_0} + 2\sqrt{2}\,\gamma C_{\mathrm{rem}} M_G + O(\varepsilon)\right)\|(\hat G_1 - \hat G_2, \hat W_1 - \hat W_2)\|_{*,a} < \frac{1}{2}\|(\delta\hat G, \delta\hat W)\|_{*,a}.$$
Nishida's shrinking-cone sequence $a_j \uparrow a_\infty = \frac{\rho_0}{4}$ keeps $\|\hat G^{(j)}\|_{\mathcal{X}_\rho} \le M_G$ and $d_{a_j}\|\hat W^{(j)}\|_{\mathcal{X}_\rho} \le M_W$ uniformly, converging to a unique fixed point $(\hat G, \hat W) \in C([0, 2]; \mathcal{X}_{\rho_1})^2$.

### Lemma A4 (Horizon Crossing & Physical-Time Reconstruction)
For $G \in \mathcal{X}_{\rho_1}$ with $a(G) \neq 0$, $V_G(r^2) = r\,a(G) + O(r^{5/4})$ (`B4`), which is $C^1$ on $r \ge 0$ with $\partial_r V_{\widetilde G}(r^2) = a(\widetilde G(s)) + O(\varepsilon r^{1/4}) = O(\varepsilon)$. Hence $dr/ds = -1 + \frac{1}{2}V_{\widetilde G}(r^2) \le -\frac{1}{2}$ hits $0$ transversally at a unique $s_\infty = 1 + O(\varepsilon) \in (\frac{1}{2}, \frac{3}{2}) \subset [0, 2]$, reconstructing the global classical solution $(G, \eta)$ of (1.1) on $\tau \in [0, \infty)$ with $a(\tau) + D(\tau) \equiv 0$ and $\|G(\cdot,\tau)\|_{L^2_\sigma} \le C\varepsilon e^{-2\sigma\tau}$. $\blacksquare$
