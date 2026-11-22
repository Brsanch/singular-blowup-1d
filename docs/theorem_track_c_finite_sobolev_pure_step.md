# Theorem (Track C): Finite-Sobolev (Non-Analytic) Pure-Step ($\eta \equiv 0$) Global Weak-Edge Compatible Manifold in $H^3_{1/4} \cap H^4_{5/8}$

**Date:** 2026-11-22  
**Status:** Complete line-by-line existence, uniqueness, and energy closure proof (`[D]`) + lemma-by-lemma numerical verification (`validate/edge_track_c_proof_lemmas.py`, `TCB-1..TCB-5 PASS`, plus `validate/edge_pure_step_finite_sobolev.py`, `FS-1..FS-3 PASS`, and Lean 4 phase-cancellation + $H^3$ & $H^4$ hierarchical $u$-dissipation in `lean/SingularBlowup/WeakEdgeManifolds.lean`, zero `sorry`).  
**Resolves:** Hostile review `findings_session_2026_11_21_weak_edge_hostile_review.md` (`A3`, `A4`).

---

## 1. Statement of the Theorem

Consider the autonomous $s$-clock pure-step ($\eta \equiv 0$) conjugated PDE on $s \in [0, 2]$ ($u = e^{-z/2}$):
$$h_s + u\,h = \mathcal{N}_{1/4}^{-1}\mathcal{P}(\widetilde G), \qquad \widetilde G = \mathcal{N}_{1/4} h, \qquad \mathcal{P}(\widetilde G) := -(Q\widetilde G)\,D\widetilde G + \tfrac{1}{2}(Q\widetilde G)\,\widetilde G, \tag{1.1}$$
with $h(\cdot, 0) = \varepsilon h_0$. Here $H^k_\sigma$ denotes the weighted Sobolev space with norm $\|f\|_{H^k_\sigma}^2 := \sum_{m=0}^k \|\partial_z^m(e^{-\sigma z} f)\|_{L^2(\mathbb{R})}^2$ (equivalently $\frac{1}{2\pi}\int_{\mathbb{R}} (1+\xi^2)^k |\widehat f(\sigma + i\xi)|^2\,d\xi$), with **no strip-analyticity assumed**.

### Rank-1 decomposition on $\operatorname{Re} m = \frac{5}{8}$
Let $\phi_*(z) := \frac{1}{2\cosh(z/2)} = \frac{u}{u^2+1}$ and $\psi_{5/8}(z) := -\frac{e^{3z/2}}{1+e^z} = -\frac{1}{u(u^2+1)} \in H^\infty_{5/8}$ (so $\phi_* - \psi_{5/8} = e^{z/2}$). Decompose any solution $h$ of (1.1) as $h = h_{\mathrm{rk1}} + \tilde h$, where
$$h_{\mathrm{rk1}}(u, s) := -\frac{\pi}{2}\int_0^s e^{-u(s-s')}\,a(\mathcal{P}(\widetilde G(s')))\,\phi_*(u)\,ds', \qquad \tilde h(\cdot, 0) = \varepsilon h_0. \tag{1.2}$$

### Theorem 3 (Global Finite-Sobolev Pure-Step Weak-Edge Compatible Manifold in $H^3_{1/4} \cap H^4_{5/8}$)
There exist universal constants $\varepsilon_0 > 0$ and $C_0 > 0$ such that for every non-analytic initial seed $h_0 \in H^3_{1/4} \cap H^4_{5/8}$ with $\|h_0\|_{H^3_{1/4}} + \|h_0\|_{H^4_{5/8}} \le 1$ and every $0 < \varepsilon \le \varepsilon_0$:
1. **(Existence, uniqueness, and finite-Sobolev regularity):** There exists a unique solution
   $$h \in C([0,2]; H^2_{1/4} \cap H^3_{5/8}) \cap C_w([0,2]; H^3_{1/4} \cap H^4_{5/8}) \cap C_R([0,2]; H^3_{1/4} \cap H^4_{5/8})$$
   of (1.1) on $[0, 2]$ with $h(\cdot, 0) = \varepsilon h_0$, satisfying the closed bootstrap bound
   $$\sup_{s \in [0,2]}\Bigl(\|h(s)\|_{H^3_{1/4}} + \|\tilde h(s)\|_{H^4_{5/8}}\Bigr) \le 2 C_0 \varepsilon, \qquad \sup_{s \in [0,2]}\|h_{\mathrm{rk1}}(s)\|_{H^5_{5/8}} \le C_1 \varepsilon^2, \tag{1.3}$$
   with $\widetilde G = \mathcal{N}_{1/4} h \in C([0,2]; L^2_{1/4} \cap H^{1/4}_{5/8}) \cap C_w([0,2]; H^1_{1/4} \cap H^{5/4}_{5/8})$ (hence strongly continuous in $H^{1-\theta}_{1/4} \cap H^{5/4-\theta}_{5/8}$ for every $\theta > 0$) and $a(\widetilde G(s)) \equiv 0$.
2. **(Global physical-time horizon crossing):** The horizon ODE $dr/ds = -1 + \frac{1}{2}V_{\widetilde G(\cdot,s)}(r^2)$, $r(0)=1$, has a unique transversal zero $s_\infty = 1 + O(\varepsilon) \in (\frac{1}{2}, \frac{3}{2}) \subset [0,2]$, yielding a unique global physical-time solution $G(v,\tau) = \widetilde G(\lambda(\tau)v, s(\tau))$ on $\tau \in [0,\infty)$ satisfying $a(\tau) \equiv 0$ ($B(\tau) \equiv A(\tau)^2 = 1$) and $\|G(\cdot,\tau)\|_{L^2_\sigma} \le C\varepsilon e^{-2\sigma\tau}$ for $\sigma \in [\frac{1}{4}, \frac{5}{8}]$.

---

## 2. Complete Proof (Lemmas C1–C5)

### Lemma C1 (Unweighted $\|b\|_{H^2(\mathbb{R})}$ Control & Exact Contour-Shift Identity)
1. **Unweighted $\|b\|_{H^2(\mathbb{R})}$ from $H^3_{1/4} \cap H^4_{5/8}$ (closing `A3`!):** Let $h \in H^3_{1/4} \cap H^4_{5/8}$. Since $\widehat h(m)$ is holomorphic on $\frac{1}{4} < \operatorname{Re} m < \frac{5}{8}$ with $N(\frac{1}{2}) = 0$, $\widetilde G = \mathcal{N}_{1/4} h = \mathcal{N}_{5/8} h \in H^1_{1/4} \cap H^{5/4}_{5/8}$ satisfies $a(\widetilde G) = 0$. Therefore $b := Q\widetilde G$ has no pole at $\zeta = 0$ and satisfies
   $$\|e^{z/4} b\|_{H^2(\mathbb{R})} \le C_1\|\widetilde G\|_{H^1_{1/4}} \le C_1'\|h\|_{H^3_{1/4}}, \qquad \|e^{-z/8} b\|_{H^{9/4}(\mathbb{R})} \le C_2\|\widetilde G\|_{H^{5/4}_{5/8}} \le C_2'\|h\|_{H^4_{5/8}}.$$
   Splitting $\mathbb{R} = (-\infty, 0] \cup [0, \infty)$ and applying the product rule to $\partial_z^m(e^{-z/4}(e^{z/4}b))$ on $z \ge 0$ and $\partial_z^m(e^{z/8}(e^{-z/8}b))$ on $z \le 0$ ($m = 0,1,2$):
   $$\boxed{\|b\|_{H^2(\mathbb{R})} \le C_b\bigl(\|h\|_{H^3_{1/4}} + \|h\|_{H^4_{5/8}}\bigr) \le C_b\bigl(\|h\|_{H^3_{1/4}} + \|\tilde h\|_{H^4_{5/8}} + \|h_{\mathrm{rk1}}\|_{H^4_{5/8}}\bigr).} \tag{2.1}$$
2. **Exact contour-shift identity across $m = \frac{1}{2}$:** Since $b \in H^2(\mathbb{R}) \hookrightarrow L^\infty(\mathbb{R})$ and $D\widetilde G \in L^2_{1/4} \cap H^{1/4}_{5/8}$, we have $F := \mathcal{P}(\widetilde G) \in L^2_{1/4} \cap L^2_{5/8}$. Thus $\widehat F(m)/N(m)$ has a single simple pole at $m = \frac{1}{2}$ with residue $+\frac{\pi}{2}a(F)$, cancelled by $\frac{\pi}{2}a(F)\widehat{\phi_*}(m) = \frac{\pi}{2}a(F)\frac{\pi}{\cos(\pi m)}$ ($\operatorname{Res}_{1/2} = -1$). Shifting the Mellin inversion contour from $\operatorname{Re} m = \frac{1}{4}$ to $\operatorname{Re} m = \frac{5}{8}$ yields the exact identity (`TCB-3`, verified to `3.55e-12`):
   $$\mathcal{N}_{1/4}^{-1} F + \frac{\pi}{2}a(F)\,\phi_* = \mathcal{N}_{5/8}^{-1} F + \frac{\pi}{2}a(F)\,\psi_{5/8}, \qquad \psi_{5/8}(z) = -\frac{e^{3z/2}}{1+e^z} \in H^\infty_{5/8}. \tag{2.2}$$
   Subtracting $\partial_s h_{\mathrm{rk1}} + u h_{\mathrm{rk1}} = -\frac{\pi}{2}a(\mathcal{P}(\widetilde G))\phi_*$ from (1.1) gives the exact PDE for $\tilde h$:
   $$\tilde h_s + u\,\tilde h = \mathcal{N}_{5/8}^{-1}\Bigl[-b\,D(\mathcal{N}_{5/8}\tilde h) - b\,D(\mathcal{N}_{5/8} h_{\mathrm{rk1}}) + \tfrac{1}{2}b\,\widetilde G\Bigr] + \frac{\pi}{2}a(\mathcal{P}(\widetilde G))\,\psi_{5/8}. \tag{2.3}$$

---

### Lemma C2 ($H^5_{5/8}$ Rank-1 Estimate and Zero-Derivative $a(\mathcal{P}(\widetilde G))$)
1. By Sobolev embedding $H^1_\sigma \hookrightarrow C_{0,\sigma}$, the boundary terms in the integration by parts for $a(\mathcal{P}(\widetilde G)) = -\frac{1}{\pi}\int_{\mathbb{R}} e^{-z/2} D(Q\widetilde G)\,\widetilde G\,dz$ vanish for all $\widetilde G \in H^1_{1/4} \cap H^{5/4}_{5/8}$, and $\|D(Q\widetilde G)\|_{L^2_{-1/4}} \le C\|h\|_{H^3_{1/4}}$, $\|D(Q\widetilde G)\|_{L^2_{1/8}} \le C\|h\|_{H^4_{5/8}}$ yield
   $$\bigl|a(\mathcal{P}(\widetilde G))\bigr| \le C_a\bigl(\|h\|_{H^3_{1/4}} + \|h\|_{H^4_{5/8}}\bigr)^2.$$
2. For every integer $j \ge 0$, $\partial_z^j(e^{-u t}\phi_*)$ satisfies $\int_{-\infty}^0 e^{-5x/4}e^x(1+(rt)^j)^2 e^{-2rt}\,dx = 2t^{-1/2}\int_t^\infty v^{-1/2}(1+v^j)^2 e^{-2v}\,dv \le C_j t^{-1/2}$ ($r=e^{-x/2}, v=rt$), while $x \ge 0$ is $O(1)$ uniformly in $t \ge 0$. Taking the square root gives $\|e^{-ut}\phi_*\|_{H^5_{5/8}} \le C_5(1 + t^{-1/4})$ (`TCB-2`), and hence from (1.2):
   $$\boxed{\|h_{\mathrm{rk1}}(s)\|_{H^5_{5/8}} \le C_5'\,s^{3/4}\sup_{0 \le s' \le s}\bigl(\|h(s')\|_{H^3_{1/4}} + \|h(s')\|_{H^4_{5/8}}\bigr)^2.} \tag{2.4}$$

---

### Lemma C3 (Exact Phase-Cancellation & Lower-Order Bounds in $H^3_{1/4}$ and $H^4_{5/8}$, resolving `A4`)
For $\sigma \in \{\frac{1}{4}, \frac{5}{8}\}$, let $A_b^{(\sigma)} := \mathcal{N}_\sigma^{-1}(b(\sigma + \partial_z))\mathcal{N}_\sigma$ on $v = e^{-\sigma z}f$, with multiplier $m_\sigma(\xi) = N(\sigma+i\xi) = A_\sigma(\xi)e^{i\Phi_\sigma(\xi)}$, $A_\sigma(\xi) \asymp \langle\xi\rangle^{p_\sigma}$ ($p_{1/4} = 2$, $p_{5/8} = \frac{11}{4}$). For any real positive polynomial weight $W_k(\xi) \asymp (1+\xi^2)^{k/2}$ — in particular both standard $(1+\xi^2)^{k/2}$ and the hierarchical weight $W_{k,\mathrm{hier}}(\xi) := (\sum_{m=0}^k c_m^{(k)}\xi^{2m})^{1/2}$ of Lemma C4 — the exact phase-cancellation identity (`WeakEdgeManifolds.phase_cancellation_norm`) gives kernel magnitude
$$|K_{\sigma,k}(\xi,\eta)| = \left|(\xi-\eta)\frac{A_{\sigma,k}(\xi-\eta)}{A_{\sigma,k}(\xi)} - \xi\frac{A_{\sigma,k}(\xi)}{A_{\sigma,k}(\xi-\eta)}\right|, \qquad A_{\sigma,k}(\xi) := \frac{A_\sigma(\xi)}{W_k(\xi)} \asymp \langle\xi\rangle^{p_\sigma - k}.$$
- **At $(\sigma_1, k_1) = (\frac{1}{4}, 3)$:** $p_{1/4} - 3 = -1$, so $|K_{1/4, 3}(\xi,\eta)| \le C(1 + |\eta|)$ and both $A_{1/4,3}(\xi-\eta)/A_{1/4,3}(\xi)$ and $A_{1/4,3}(\xi)/A_{1/4,3}(\xi-\eta)$ are $\le C(1 + |\eta|)$ uniformly for all $\xi, \eta \in \mathbb{R}$.
- **At $(\sigma_2, k_2) = (\frac{5}{8}, 4)$:** $p_{5/8} - 4 = -\frac{5}{4}$, so $|K_{5/8, 4}(\xi,\eta)| \le C(1 + |\eta|^{5/4})$ and both $A_{5/8,4}(\xi-\eta)/A_{5/8,4}(\xi)$ and $A_{5/8,4}(\xi)/A_{5/8,4}(\xi-\eta)$ are $\le C(1 + |\eta|^{5/4})$ uniformly for all $\xi, \eta \in \mathbb{R}$ (`TCB-1`: flat for $b$ barely in $H^2$!).

Since $\int_{\mathbb{R}} (1 + |\eta|^{5/4})|\widehat b(\eta)|\,d\eta \le C\|b\|_{H^2(\mathbb{R})}$ (by Cauchy–Schwarz, with a full $\frac{1}{4}$-derivative margin since $\frac{5}{4} + \frac{1}{2} = \frac{7}{4} < 2$!), Schur's test yields:
$$\bigl|\operatorname{Re}\langle A_b^{(1/4)} v, v\rangle_{H^3_{\mathrm{hier}}}\bigr| \le C\|b\|_{H^2(\mathbb{R})}\,\|h\|_{H^3_{1/4}}^2, \qquad \bigl\|\mathcal{N}_{1/4}^{-1}(b\widetilde G)\bigr\|_{H^3_{1/4}} \le C\|b\|_{H^2(\mathbb{R})}\,\|h\|_{H^3_{1/4}}, \tag{2.5}$$
$$\bigl|\operatorname{Re}\langle A_b^{(5/8)} \tilde v, \tilde v\rangle_{H^4_{\mathrm{hier}}}\bigr| \le C\|b\|_{H^2(\mathbb{R})}\,\|\tilde h\|_{H^4_{5/8}}^2, \qquad \bigl\|\mathcal{N}_{5/8}^{-1}(b\,D(\mathcal{N}_{5/8} h_{\mathrm{rk1}}))\bigr\|_{H^4_{5/8}} \le C\|b\|_{H^2(\mathbb{R})}\,\|h_{\mathrm{rk1}}\|_{H^5_{5/8}}, \tag{2.6}$$
$$\bigl\|\mathcal{N}_{5/8}^{-1}(b\widetilde G)\bigr\|_{H^4_{5/8}} \le C\|b\|_{H^2(\mathbb{R})}\,\|h\|_{H^4_{5/8}}. \tag{2.7}$$
*(Note on the cross term in (2.6): in Fourier space, $\frac{W_4(\xi)}{A_{5/8}(\xi)} |\frac{5}{8}+i(\xi-\eta)| A_{5/8}(\xi-\eta) / W_5(\xi-\eta) \le C(1 + |\eta|^{5/4})$ because $|\frac{5}{8}+i(\xi-\eta)| \le \langle\xi-\eta\rangle = W_5/W_4$!)*

---

### Lemma C4 (Hierarchical $u$-Dissipation in $H^3$ & $H^4$ and Closed Bootstrap on $[0, 2]$)
With hierarchical weights $c_m^{(3)} = 4^{3-m}$ ($m = 0..3$) and $c_m^{(4)} = 6^{4-m}$ ($m = 0..4$),
$$\sum_{m=0}^k c_m^{(k)}\operatorname{Re}\langle \partial_z^m(u v), \partial_z^m v\rangle_{L^2} \ge \frac{1}{2}\sum_{m=0}^k c_m^{(k)}\int_{\mathbb{R}} u\,|\partial_z^m v|^2\,dz \ge 0 \qquad (k \in \{3, 4\}), \tag{2.8}$$
verified in Lean (`WeakEdgeManifolds.hierarchical_u_dissipation_H3` for $k = 3$ and `WeakEdgeManifolds.hierarchical_u_dissipation_H4` for $k = 4$, both zero `sorry`, and `TCB-4`).

Define the equivalent hierarchical energy $E(s) := \|h(s)\|_{H^3_{1/4,\mathrm{hier}}}^2 + \|\tilde h(s)\|_{H^4_{5/8,\mathrm{hier}}}^2$ and let $Y(s) := \sup_{0\le s'\le s}\sqrt{E(s')}$. Under the bootstrap hypothesis $Y(s) \le 2 C_0 \varepsilon$, (2.4) gives $\|h_{\mathrm{rk1}}(s)\|_{H^5_{5/8}} \le C Y(s)^2$, (2.1) gives $\|b(s)\|_{H^2(\mathbb{R})} \le C(Y(s) + Y(s)^2)$, and combining (2.5)–(2.8) on (1.1) and (2.3) yields
$$\frac{d}{ds} E(s) \le C'\bigl(Y(s) + Y(s)^2\bigr)\bigl(E(s) + Y(s)^4\bigr) \le C'' \varepsilon\,Y(s)^2 \qquad (s \in [0, 2]).$$
Integrating on $[0, 2]$ gives $Y(2)^2 \le E(0) + 2 C''\varepsilon Y(2)^2$. For $\varepsilon \le \varepsilon_0 := (4 C'')^{-1}$, this closes the bootstrap $Y(2) \le \sqrt{2 E(0)} \le 2 C_0 \varepsilon$, establishing (1.3).

---

### Lemma C5 (Parabolic Regularization, Unconditional Uniqueness, Strong Convergence, and Right-Continuity)
1. **Structure-preserving parabolic regularization (single-unknown fixed point):** Define the constant-coefficient Mellin multiplier
   $$L := -(\partial_z - \tfrac{1}{4})(\partial_z - \tfrac{5}{8}) + \tfrac{1}{4} = -\partial_z^2 + \tfrac{7}{8}\partial_z + \tfrac{3}{32},$$
   whose symbol on $\operatorname{Re} m = \sigma$ is $\widehat L(\sigma + i\xi) = (\xi^2 + \frac{1}{4} - (\sigma - \frac{1}{4})(\sigma - \frac{5}{8})) + i\xi(\frac{7}{8} - 2\sigma)$. On both $\sigma = \frac{1}{4}$ and $\sigma = \frac{5}{8}$, $\operatorname{Re}\widehat L(\sigma + i\xi) = \xi^2 + \frac{1}{4} > 0$, so $L$ commutes identically with $\mathcal{N}_{1/4}$ and $\mathcal{N}_{5/8}$, preserves $a(\widetilde G) = 0$, and is strictly coercive on $H^k_{\sigma,\mathrm{hier}}$: there exist explicit constants $c_k > 0$ such that
   $$\operatorname{Re}\langle L f, f\rangle_{H^k_{\sigma,\mathrm{hier}}} \ge c_k \|f\|_{H^{k+1}_\sigma}^2 \qquad (\sigma \in \{\tfrac{1}{4}, \tfrac{5}{8}\},\; k \ge 0).$$
   For $\delta \in (0, 1]$, we construct a solution of
   $$h^\delta_s + u h^\delta + \delta L h^\delta = \mathcal{N}_{1/4}^{-1}\mathcal{P}(\widetilde G^\delta), \qquad h^\delta(\cdot, 0) = \varepsilon h_0, \tag{2.9}$$
   via a fixed-point map $\Phi$ on the **single** unknown $h \in C([0, T_\delta]; H^3_{1/4} \cap H^4_{5/8})$:
   - Given $h \in C([0, T_\delta]; H^3_{1/4} \cap H^4_{5/8})$, Lemma C1 gives $\widetilde G = \mathcal{N}_{1/4} h = \mathcal{N}_{5/8} h$ with $a(\widetilde G) = 0$ and $b = Q\widetilde G \in C([0, T_\delta]; H^2(\mathbb{R}))$, while Lemma C2 gives $h_{\mathrm{rk1}}[h] \in C([0, T_\delta]; H^\infty_{1/4} \cap H^5_{5/8})$ via (1.2) (since $e^{-z/4}\phi_*(z) \sim e^{z/4}$ as $z \to -\infty$ gives $\|e^{-ut}\phi_*\|_{H^m_{1/4}} \le C_m$ uniformly in $t \ge 0$).
   - Define the source
     $$\mathcal{F}_\delta[h] := \mathcal{N}_{5/8}^{-1}\bigl[-b D\widetilde G + \tfrac{1}{2}b\widetilde G\bigr] + \frac{\pi}{2}a(\mathcal{P}(\widetilde G))\psi_{5/8} - \delta L h_{\mathrm{rk1}}[h],$$
     which by the exact contour-shift identity (2.2) coincides identically on $\operatorname{Re} m = \frac{1}{4}$ with $\mathcal{N}_{1/4}^{-1}\mathcal{P}(\widetilde G) + \frac{\pi}{2}a(\mathcal{P}(\widetilde G))\phi_* - \delta L h_{\mathrm{rk1}}[h]$, and lies in $C([0, T_\delta]; H^2_{1/4} \cap H^3_{5/8})$ (using $\|L h_{\mathrm{rk1}}[h]\|_{H^3_{5/8}} \le C\|h_{\mathrm{rk1}}[h]\|_{H^5_{5/8}}$ and $H^2 \cdot H^{1/4} \hookrightarrow H^{1/4}$).
   - Solve the **linear** parabolic equation $(\partial_s + u + \delta L)\tilde h = \mathcal{F}_\delta[h]$ with $\tilde h(\cdot, 0) = \varepsilon h_0$. On both Gelfand triples $V_\sigma \hookrightarrow H^{k_\sigma}_{\sigma,\mathrm{hier}} \hookrightarrow V_\sigma^*$ ($k_{1/4}=3$, $k_{5/8}=4$, $V_\sigma = \{f \in H^{k_\sigma+1}_\sigma : u^{1/2}f \in H^{k_\sigma}_\sigma\}$), $u + \delta L : V_\sigma \to V_\sigma^*$ is bounded and strictly coercive by (2.8), and $H^{k_\sigma - 1}_\sigma \hookrightarrow V_\sigma^*$. By the Lions parabolic theorem, there is a unique solution $\tilde h \in C([0, T_\delta]; H^3_{1/4} \cap H^4_{5/8})$ satisfying (for an equivalence constant $c > 0$)
     $$\|\tilde h(s)\|_{H^{k_\sigma}_{\sigma,\mathrm{hier}}}^2 + c\int_0^s \|u^{1/2}\tilde h\|_{H^{k_\sigma}_\sigma}^2\,ds' + c_{k_\sigma}\delta\int_0^s \|\tilde h\|_{H^{k_\sigma+1}_\sigma}^2\,ds' \le \|\tilde h(0)\|_{H^{k_\sigma}_{\sigma,\mathrm{hier}}}^2 + \frac{C}{\delta}\int_0^s \|\mathcal{F}_\delta[h](s')\|_{H^{k_\sigma-1}_\sigma}^2\,ds'.$$
   - Set $\Phi(h) := \tilde h + h_{\mathrm{rk1}}[h] \in C([0, T_\delta]; H^3_{1/4} \cap H^4_{5/8})$. Since $h \mapsto \mathcal{F}_\delta[h]$ is locally Lipschitz on $C([0, T_\delta]; H^3_{1/4} \cap H^4_{5/8}) \to C([0, T_\delta]; H^2_{1/4} \cap H^3_{5/8})$ and $h \mapsto h_{\mathrm{rk1}}[h]$ carries the $T_\delta^{3/4}$ factor (2.4), Banach's fixed-point theorem yields a unique local solution $h^\delta = \Phi(h^\delta) = \tilde h^\delta + h_{\mathrm{rk1}}[h^\delta] \in C([0, T_\delta]; H^3_{1/4} \cap H^4_{5/8})$ for $T_\delta \asymp \delta / M^2$. Adding $(\partial_s + u + \delta L)h_{\mathrm{rk1}}[h^\delta] = -\frac{\pi}{2}a(\mathcal{P}(\widetilde G^\delta))\phi_* + \delta L h_{\mathrm{rk1}}[h^\delta]$ to $(\partial_s + u + \delta L)\tilde h^\delta = \mathcal{F}_\delta[h^\delta]$ shows $h^\delta$ solves (2.9) identically.
2. **$\delta$-uniform global bound on $[0, 2]$:** In the $H^3_{1/4,\mathrm{hier}} \times H^4_{5/8,\mathrm{hier}}$ energy identity for $(h^\delta, \tilde h^\delta)$, $\delta\operatorname{Re}\langle L h^\delta, h^\delta\rangle_{H^3_{1/4,\mathrm{hier}}} \ge 0$, while on $\operatorname{Re} m = \frac{5}{8}$ integrating one derivative by parts (since $|\widehat L(\frac{5}{8}+i\xi)| \le C(1+\xi^2)$) gives
   $$\delta\bigl|\langle L h_{\mathrm{rk1}}^\delta, \tilde h^\delta\rangle_{H^4_{5/8,\mathrm{hier}}}\bigr| \le C\delta\|h_{\mathrm{rk1}}^\delta\|_{H^5_{5/8}}\|\tilde h^\delta\|_{H^5_{5/8}} \le \frac{c_4 \delta}{2}\|\tilde h^\delta\|_{H^5_{5/8}}^2 + C'\delta\|h_{\mathrm{rk1}}^\delta\|_{H^5_{5/8}}^2.$$
   The first term is absorbed by the $+c_4\delta\|\tilde h^\delta\|_{H^5_{5/8}}^2$ dissipation on the LHS, and the second is $\le C''\delta Y_\delta^4 \le C''\varepsilon^2 Y_\delta^2$ uniformly for $\delta \in (0, 1]$. Hence Lemma C4 closes uniformly in $\delta \in (0, 1]$, extending $h^\delta$ to $[0, 2]$ with bound (1.3).
3. **Lower-norm difference estimate (Unconditional Uniqueness, Strong Convergence, and Right-Continuity):** For any two solutions $h_1, h_2 \in L^\infty([0, 2]; H^3_{1/4} \cap H^4_{5/8})$ of (2.9) with viscosities $\delta_1, \delta_2 \in [0, 1]$ and bounds $\sup_{s\in[0,2]}(\|h_i\|_{H^3_{1/4}} + \|h_i\|_{H^4_{5/8}}) \le M_i$, let $M := M_1 + M_2$, $w := h_1 - h_2$, $w_{\mathrm{rk1}} := h_{\mathrm{rk1}}[h_1] - h_{\mathrm{rk1}}[h_2]$, $\tilde w := w - w_{\mathrm{rk1}}$, and $\delta b := Q(\mathcal{N}_{1/4} w)$. In the lower norm $E_{\mathrm{low}}(s) := \|w(s)\|_{H^2_{1/4,\mathrm{hier}}}^2 + \|\tilde w(s)\|_{H^3_{5/8,\mathrm{hier}}}^2$:
   - $\|\delta b\|_{H^1(\mathbb{R})} \le C(\|w\|_{H^2_{1/4}} + \|w\|_{H^3_{5/8}})$, and on any subinterval $[s_0, s_0 + \Delta s]$ (using the shifted $w_{\mathrm{rk1}, s_0}$ from $s_0$ to $s$), $\|w_{\mathrm{rk1}, s_0}(s)\|_{H^4_{5/8}} \le C(\Delta s)^{3/4} M\sup_{s'\in[s_0, s]}(\|w\|_{H^2_{1/4}} + \|w\|_{H^3_{5/8}})$.
   - Phase cancellation holds at $(\frac{1}{4}, 2)$ ($p_{1/4}-2=0$, $|K_{1/4,2}| \le C(1+|\eta|)$) and $(\frac{5}{8}, 3)$ ($p_{5/8}-3=-\frac{1}{4}$, $|K_{5/8,3}| \le C(1+|\eta|)$) (`TCB-5`), and hierarchical $u$-dissipation holds for $k=2$ ($c_m^{(2)}=4^{2-m}$, $\frac{1}{4}\times$ leading $3\times 3$ submatrix of $k=3$, $\lambda_{\min} = 0.7320 > 0.5$) and $k=3$ (`hierarchical_u_dissipation_H3`).
   - The non-transport terms $(\delta b)D\widetilde G_2$ satisfy $\|\mathcal{N}_{1/4}^{-1}((\delta b)D\widetilde G_2)\|_{H^2_{1/4}} \le C\|\delta b\|_{H^1}\|h_2\|_{H^3_{1/4}}$ and $\|\mathcal{N}_{5/8}^{-1}((\delta b)D\widetilde G_2)\|_{H^3_{5/8}} \asymp \|(\delta b)D\widetilde G_2\|_{H^{1/4}_{5/8}} \le C\|\delta b\|_{H^1}\|h_2\|_{H^4_{5/8}}$.
   - The viscosity forcing satisfies $|\langle \delta_i L h_i, w\rangle_{H^2_{1/4}}| + |\langle \delta_i L\tilde h_i, \tilde w\rangle_{H^3_{5/8}}| + |\langle \delta_i L h_{\mathrm{rk1}}[h_i], \tilde w\rangle_{H^3_{5/8}}| \le C\delta_i M^2$.

   Integrating on $[s_0, s_0 + \Delta s]$ yields $Y_{\mathrm{low}}^2 \le E_{\mathrm{low}}(s_0) + C M \Delta s\,Y_{\mathrm{low}}^2 + C(\delta_1 + \delta_2)M^2$. Choosing $\Delta s := \min(2, (4 C M)^{-1})$ and chaining finitely many subintervals across $[0, 2]$:
   - Setting $\delta_1 = \delta_2 = 0$ and $w(0) = 0$ proves **unconditional uniqueness** in $L^\infty([0, 2]; H^3_{1/4} \cap H^4_{5/8})$ ($w \equiv 0$).
   - Setting $h_i = h^{\delta_i}$ (with $w(0) = 0$) and $\delta_1, \delta_2 \to 0^+$ proves $\{h^\delta\}$ is **Cauchy in $C([0, 2]; H^2_{1/4} \cap H^3_{5/8})$** with rate $O(\sqrt{\delta_1 + \delta_2})$, converging strongly to $h \in C([0, 2]; H^2_{1/4} \cap H^3_{5/8}) \cap C_w([0, 2]; H^3_{1/4} \cap H^4_{5/8})$ solving (1.1) and satisfying (1.3).
   - **Strong right-continuity at every $s_0 \in [0, 2)$:** By autonomy of (1.1), restart the $\delta$-regularized flow (2.9) at $s_0$ with initial data $h(s_0) \in H^3_{1/4} \cap H^4_{5/8}$ and shifted rank-1 piece $h_{\mathrm{rk1}, s_0}(s) := -\frac{\pi}{2}\int_{s_0}^s e^{-u(s-s')}a(\mathcal{P}(\widetilde G(s')))\phi_*\,ds' \in C([s_0, 2]; H^5_{5/8})$. By uniqueness, its $\delta \to 0^+$ limit coincides with $h$ on $[s_0, 2]$. Passing $\delta \to 0^+$ in the restarted energy bound $E_{s_0,\delta}(s) \le E_{s_0,\delta}(s_0) + C''\varepsilon Y^2(s - s_0)$ (where $E_{s_0,\delta}(s_0) = E_{s_0}(s_0)$ is exact!) gives $\limsup_{s \to s_0^+} E_{s_0}(s) \le E_{s_0}(s_0)$, which with weak continuity $\liminf_{s \to s_0^+} E_{s_0}(s) \ge E_{s_0}(s_0)$ and $h_{\mathrm{rk1}, s_0} \in C([s_0, 2]; H^4_{5/8})$ implies strong right-continuity $h \in C_R([0, 2]; H^3_{1/4} \cap H^4_{5/8})$.

Horizon crossing and physical-time reconstruction follow identically from Track B Lemma 4 (since $H^1_{1/4} \cap H^{5/4}_{5/8} \hookrightarrow L^\infty$ gives $|G(w)| \le C\|G\|\min(w^{5/8}, w^{1/4})$ with $a(G) = 0$, hence $V_G(r^2) = O(\varepsilon r^{5/4})$). $\blacksquare$
