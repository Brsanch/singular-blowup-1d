# Finite-Sobolev (Non-Analytic) Energy Closure for the Pure-Step ($\eta \equiv 0$) Global Weak-Edge Manifold (2026-11-21)

**Context (updated 2026-11-22 after hostile review `findings_session_2026_11_21_weak_edge_hostile_review.md`):** Investigates whether the pure-step ($\eta \equiv 0$) conjugated equation
$$h_s + u\,h = \mathcal{N}_g^{-1}\Bigl[-(Q\widetilde G)\,D\widetilde G + \tfrac{1}{2}(Q\widetilde G)\,\widetilde G\Bigr], \qquad \widetilde G = \mathcal{N}_g h, \quad u = e^{-z/2},$$
closes in finite Sobolev space without strip-analyticity. **Hostile-review verdict (A3, A4):**
- **What holds (`[D, L, N]`):** The exact phase-cancellation identity (Lean-verified in `WeakEdgeManifolds.phase_cancellation_norm`), the top-order symmetric-part bound $|\operatorname{Re}\langle A_b v, v\rangle_{H^k}| \le C_k\|b\|_{H^2}\|v\|_{H^k}^2$ for **$k \in \{2, 3\}$ only** at $\sigma = \frac{1}{4}$, the hierarchical $H^3$ $u$-dissipation coercivity (Lean-verified in `hierarchical_u_dissipation_H3`), and filter-independent monotone-decreasing $H^2, H^3$ evolution of a non-analytic $C^\infty_c$ bump on clean domains `[-72, 72]`.
- **What does NOT close as stated (A3, A4):** (1) The $\|b\|_{H^2}$ bound is **false at $k = 0$** (requires $\|b\|_{H^4}$) and unproven at $k = 1$ (requires $\|b\|_{H^3}$); (2) $h \in L^\infty_s H^3_{1/4}$ plus $u$-dissipation $L^2_s H^3_{1/2}$ falls short by $\delta > 0$ derivatives at the edge of controlling unweighted $\|b\|_{H^{3/2+\delta}}$ (since $b = -2\mathcal{N}(uh) \in H^j(\mathbb{R}) \iff h \in H^{j+3/2}_{1/2}$). Closing the second-weight $\sigma_2 > \frac{1}{2}$ rank-1 defect estimate remains **open**.

**Evidence tiers:** `[D, L]` for the phase-cancellation identity, $k \in \{2,3\}$ symmetric-part bound, and $H^3$ $u$-dissipation; `[C]` / **OPEN** for full finite-Sobolev nonlinear closure (A3); `[N]` for `validate/edge_pure_step_finite_sobolev.py` (`FS-1..FS-3 PASS`).

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
Hence, for **$k \in \{2, 3\}$** at $\sigma = \frac{1}{4}$ ($p = 2$, where $p - k \in \{0, -1\}$; note that $k = 0$ requires $\|b\|_{H^4}$ and $k = 1$ requires $\|b\|_{H^3}$, per hostile review A4):
$$\left|\operatorname{Re}\langle A_b v, v \rangle_{H^k}\right| \le C_k\,\|\widehat{b_z}\|_{L^1}\,\|v\|_{H^k}^2 \le C_k'\,\|b\|_{H^2}\,\|v\|_{H^k}^2 \qquad \textbf{(ZERO derivative loss at } k \in \{2, 3\}\textbf{)}.$$

**Check (`FS-1`):** Algebraic identity holds to `< 1e-12` (Lean-verified in `WeakEdgeManifolds.phase_cancellation_norm`); operator norm $\|(A_b + A_b^*)/2\|_{H^k \to H^k}$ is $N$-independent for $k \in \{2, 3\}$ even for $b$ barely in $H^2$ (while $k = 0$ grows unboundedly with cutoff, confirming A4).

---

## 2. Hierarchical $u$-Dissipation and Finite-Sobolev Closure (`[D, N]`, `FS-2`, `FS-3`)

Because $\partial_z u = -\frac{1}{2} u$, differentiating $u v$ yields $\partial_z^m(u v) = u(\partial_z - \frac{1}{2})^m v$, where every lower-order commutator carries the same positive weight $u = e^{-z/2}$. With hierarchical weights $c_m = 4^{3-m}$ ($m = 0,1,2,3$), Young's inequality gives
$$\sum_{m=0}^3 c_m \operatorname{Re}\bigl\langle \partial_z^m(u v),\, \partial_z^m v \bigr\rangle_{L^2} \ge \frac{1}{2}\sum_{m=0}^3 c_m \int_{\mathbb{R}} u\,|\partial_z^m v|^2\,dz \ge 0$$
(`FS-2`: min ratio `0.903016 > 0.5` on test functions; Lean-verified in `WeakEdgeManifolds.hierarchical_u_dissipation_H3`). **Note (2026-11-22):** The `A3` edge-regularity gap noted here is now closed in `docs/theorem_track_c_finite_sobolev_pure_step.md` (`H^3_{1/4} \cap H^4_{5/8}`, `TCB-1..5 PASS`, plus Lean 4 `hierarchical_u_dissipation_H4`).

**Check (`FS-3`):** Evolving a **non-analytic $C^\infty_c$ bump** $h_0(u)$ supported on $u \in (0.5, 4.0)$ through the full nonlinear pure-step ($\eta \equiv 0$) PDE on $s \in [0, 1.2]$ ($> s_\infty$) via sponge-stabilized ETD-RK2 across `([-64, 64], 4096/8192)` and `([-80, 80], 4096)` gives converged norms across all four orders (`L^2 = 1.703e-3`, `H^1 = 2.935e-3`, `H^2 = 6.986e-3`, `H^3 = 2.577e-2`, `Nz_drift = 4.3e-7`, `dom_drift = 7.7e-4`) and tail-corrected real-space quadrature `|a_quad| = 1.66e-10`.
