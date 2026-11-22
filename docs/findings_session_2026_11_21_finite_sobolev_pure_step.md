# Finite-Sobolev (Non-Analytic) Energy Closure for the Pure-Step ($\eta \equiv 0$) Global Weak-Edge Manifold (2026-11-21)

**Context:** Pushes [findings_session_2026_11_21_pure_step_nonlinear_manifold.md](findings_session_2026_11_21_pure_step_nonlinear_manifold.md) beyond strip-analyticity. Shows that for the pure-step ($\eta \equiv 0$) conjugated equation
$$h_s + u\,h = \mathcal{N}_g^{-1}\Bigl[-(Q\widetilde G)\,D\widetilde G + \tfrac{1}{2}(Q\widetilde G)\,\widetilde G\Bigr], \qquad \widetilde G = \mathcal{N}_g h, \quad u = e^{-z/2},$$
**strip-analyticity is NOT needed**: an exact **phase-cancellation identity** eliminates the logarithmic phase derivative of $\mathcal{N}_g$ completely in the symmetric part of the conjugated transport operator, yielding a **zero-derivative-loss finite-Sobolev $H^k_\sigma$ energy estimate** ($k = 3$ at $\sigma = \frac{1}{4}$, i.e. $\widetilde G \in H^1_{1/4}$).

**Evidence tiers:** `[D]` hand derivation; `[N]` numerical verification (`validate/edge_pure_step_finite_sobolev.py`, `FS-1..FS-3 PASS`). Not Lean-formalized.

---

## 1. Exact Phase-Cancellation Identity (`[D, N]`, `FS-1`)

Linearizing around a real velocity $b(z) := (Q\widetilde G)(e^z)$, the principal (highest-derivative) conjugated operator on $v := e^{-\sigma z} h \in L^2(\mathbb{R})$ is
$$A_b := \mathcal{N}_g^{-1} \circ \bigl(b(z)(\sigma + \partial_z)\bigr) \circ \mathcal{N}_g,$$
where $\mathcal{N}_g$ has Fourier multiplier $m(\xi) := N(\sigma + i\xi) = A(\xi)\,e^{i\Phi(\xi)}$ with
$$A(\xi) = |N(\sigma + i\xi)| \sim C\,|\xi|^{2\sigma + 3/2}, \qquad \Phi(\xi) = 2\xi\ln(2|\xi|) - 2\xi + O(1).$$
Although $\mathcal{N}_g \notin S^p_{1,0}$ (since $\Phi'(\xi) \sim 2\ln(2|\xi|) \to \infty$), in the symmetric part $\frac{1}{2}(A_b + A_b^*)$ on $H^k$ (weight $W_k(\xi) = (1+\xi^2)^{k/2}$, real and positive, so $m_k(\xi) := m(\xi)/W_k(\xi) = A_k(\xi) e^{i\Phi(\xi)}$ has the **identical phase** $\Phi(\xi)$!), the kernel between frequencies $\xi$ and $\xi - \eta$ is
$$K_k(\xi, \eta) := i(\xi - \eta)\frac{m_k(\xi - \eta)}{m_k(\xi)} - i\xi\frac{\overline{m_k(\xi)}}{\overline{m_k(\xi - \eta)}} = i\,e^{i(\Phi(\xi - \eta) - \Phi(\xi))}\left[(\xi - \eta)\frac{A_k(\xi - \eta)}{A_k(\xi)} - \xi\frac{A_k(\xi)}{A_k(\xi - \eta)}\right].$$
**The oscillatory phase $e^{i(\Phi(\xi-\eta)-\Phi(\xi))}$ factors out identically with modulus $1$**, and the leading $O(\xi)$ terms cancel:
$$|K_k(\xi, \eta)| = \left|(\xi - \eta)\frac{A_k(\xi - \eta)}{A_k(\xi)} - \xi\frac{A_k(\xi)}{A_k(\xi - \eta)}\right| \le C_k\,|\eta| \qquad (k \in \{2, 3\},\; \sigma = \tfrac{1}{4},\; p = 2\sigma + \tfrac{3}{2} = 2).$$
Hence, for $k \in \{0, 1, 2, 3\}$ at $\sigma = \frac{1}{4}$:
$$\left|\operatorname{Re}\langle A_b v, v \rangle_{H^k}\right| \le C_k\,\|\widehat{b_z}\|_{L^1}\,\|v\|_{H^k}^2 \le C_k'\,\|b\|_{H^2}\,\|v\|_{H^k}^2 \qquad \textbf{(ZERO derivative loss!)}.$$

**Check (`FS-1`):** Algebraic identity holds to `< 1e-12`; operator norm $\|(A_b + A_b^*)/2\|_{H^k \to H^k}$ is $N$-independent to `< 1e-4` from $N=256$ to $N=512$ (`H^0: 3.6546`, `H^1: 1.9331`, `H^2: 1.0331`, `H^3: 0.6026`) while $\|A_b\|_2$ doubles linearly with $\xi_{\max}$.

---

## 2. Hierarchical $u$-Dissipation and Finite-Sobolev Closure (`[D, N]`, `FS-2`, `FS-3`)

Because $\partial_z u = -\frac{1}{2} u$, differentiating $u v$ yields $\partial_z^m(u v) = u(\partial_z - \frac{1}{2})^m v$, where every lower-order commutator carries the same positive weight $u = e^{-z/2}$. With hierarchical weights $c_m = 4^{3-m}$ ($m = 0,1,2,3$), Young's inequality gives
$$\sum_{m=0}^3 c_m \operatorname{Re}\bigl\langle \partial_z^m(u v),\, \partial_z^m v \bigr\rangle_{L^2} \ge \frac{1}{2}\sum_{m=0}^3 c_m \int_{\mathbb{R}} u\,|\partial_z^m v|^2\,dz \ge 0$$
(`FS-2`: min ratio `> 0.5`). Combined with the exact phase-cancellation bound at $k = 3$ ($h \in H^3_{1/4} \iff \widetilde G = \mathcal{N}_g h \in H^1_{1/4} \implies b = Q\widetilde G \in H^2$ with $a(\widetilde G) \equiv 0$), the nonlinear pure-step energy closes on $s \in [0, s_\infty)$ in finite Sobolev regularity **without strip-analyticity** (with the $z \le 0$ bound on $b$ supplied by the $u$-dissipation / $H^4_{3/4}$ control after splitting the single rank-1 residue $f_{\mathcal{P}} = -\frac{\pi}{2}a(\mathcal{P})\frac{u}{u^2+1} + \tilde f_{\mathcal{P}}$).

**Check (`FS-3`):** Evolving a **non-analytic $C^\infty_c$ bump** $h_0(u)$ supported on $u \in (0.5, 4.0)$ through the full nonlinear pure-step ($\eta \equiv 0$) PDE on $s \in [0, 1.2]$ ($> s_\infty$) via ETD-RK2 at $N_z \in \{1024, 2048, 4096\}$ gives grid-independent norms ($\|h(1.2)\|_{L^2} = 1.7027\times 10^{-3}$, $\|h(1.2)\|_{H^1} = 2.9353\times 10^{-3}$, drift `< 2e-4`) and global weak-edge compatibility $\max_s |a(G(s))| = 1.76\times 10^{-6}$ with zero $\Theta$-cusp.
