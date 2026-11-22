# Adversarial Review Packet — Global Nonlinear Weak-Edge Compatibility Manifolds (2026-11-21)

**Purpose:** Self-contained review target for an independent adversarial review (e.g. Fabel) of the 2026-11-21 global nonlinear weak-edge compatibility claims in `/Volumes/4TB SD/ClaudeCode/singular-blowup-1d`.

**Git snapshot:** `main` commits `523c6a2`, `e44437a`, `dd04bb6` (built on `b4be3dc`, merged alongside `a73dcbc`).

---

## 1. Files & Reproducible Gates Under Review

| Track | Derivation Note | Verification Gate | Repro Command (via `reslock`) |
|---|---|---|---|
| **Track A:** $\sqrt{v}$ $\Theta$-cusp strip-analytic manifold | `docs/findings_session_2026_11_21_nonlinear_cusp_manifold.md` | `validate/edge_nonlinear_cusp_manifold.py` (`NM-1..NM-3`) | `/Users/bryan/.claude/bin/reslock run -- python3 validate/edge_nonlinear_cusp_manifold.py` |
| **Track B:** Pure-step ($\eta \equiv 0$) strip-analytic manifold | `docs/findings_session_2026_11_21_pure_step_nonlinear_manifold.md` | `validate/edge_pure_step_nonlinear_manifold.py` (`PM-1..PM-3`) | `/Users/bryan/.claude/bin/reslock run -- python3 validate/edge_pure_step_nonlinear_manifold.py` |
| **Track C:** Pure-step ($\eta \equiv 0$) finite-Sobolev ($H^3_{1/4}$) closure | `docs/findings_session_2026_11_21_finite_sobolev_pure_step.md` | `validate/edge_pure_step_finite_sobolev.py` (`FS-1..FS-3`) | `/Users/bryan/.claude/bin/reslock run -- python3 validate/edge_pure_step_finite_sobolev.py` |

**Background / context files in repo:**
- `paper/stability_spectra_singular_houluo.tex` (`eq:dyn`, `eq:comoving_exact`, `rem:weak_edge_jump`)
- `docs/findings_session_2026_11_21_sqrt_theta_cusp.md` + `validate/edge_sqrt_theta_cusp_closure.py` (`TC-1..3`)
- `docs/findings_session_2026_11_21_pure_step_quadratic_closure.md` + `validate/edge_pure_step_quadratic_closure.py` (`PQ-1..3`)

---

## 2. Exact Claims to Audit (with Evidence Tiers)

Every claim below is tagged `[D]` (hand derivation, not Lean-formalized) and/or `[N]` (numerical gate).

### Claim 1 — Exact Autonomous Finite-Horizon $s$-Clock Reduction (`[D, N]`, `NM-1`, `NM-3`, `PM-3`)
In the gauge $A \equiv B \equiv \xi \equiv 1$, $c_l(\tau) = 2 - V_G(\tau)$, $c_\omega(\tau) = -c_l(\tau)/2$ ($q = c_l + 2c_\omega \equiv 0$), substituting $\Omega(1+v,\tau) = v^{-1/2}(1+G(v,\tau))\mathbf{1}_{v>0}$ and $\Theta(1+v,\tau) = (\frac{\pi}{2}+\eta(v,\tau))\mathbf{1}_{v>0}$ into `(eq:dyn)` gives on $v > 0$:
$$G_\tau + (2 - V_G + QG)\,v G_v - \tfrac{1}{2}(QG)(1+G) = \sqrt{v}\,\eta_v, \qquad \eta_\tau + (2 - V_G + QG)\,v \eta_v = 0. \tag{$\star$}$$
Setting $\lambda(\tau) := \exp(-\int_0^\tau(2 - V_G)\,d\tau')$, $s(\tau) := \int_0^\tau\sqrt{\lambda(\tau')}\,d\tau' \in [0, s_\infty)$, $w := \lambda(\tau)v$, $\widetilde G(w,s) := G(v,\tau)$, $\sqrt{w}\,\widetilde H(w,s) := \eta(v,\tau)$ eliminates $V_G(\tau)$ identically by dilation covariance of $Q$:
$$\widetilde G_s + (Q\widetilde G)\,D\widetilde G = \tfrac{1}{2}(Q\widetilde G)(1+\widetilde G) + (\tfrac{1}{2}+D)\widetilde H, \qquad \widetilde H_s + (Q\widetilde G)(\tfrac{1}{2}+D)\widetilde H = 0, \qquad D := w\partial_w. \tag{$\dagger$}$$
A posteriori, $r(s) := \sqrt{\lambda(s)}$ solves $dr/ds = -1 + \frac{1}{2}V_{\widetilde G(\cdot,s)}(r^2)$, $r(0) = 1$, hitting $0$ at $s_\infty = 1 + O(\varepsilon)$, with $\tau(s) = \int_0^s ds'/r(s') \to +\infty$ as $s \to s_\infty^-$.

*(Sub-claim 1b: corrects two inaccuracies in `findings_session_2026_11_21_sqrt_theta_cusp.md` §2: quadratic $G$-source is $+\frac{1}{2}(QG_0)G_0$ not $+\frac{1}{2}(QG_0-a_0)G_0$, and $V_{G_1}(\tau) \neq e^{-\tau} V_{G_0}$.)*

### Claim 2 — $\sqrt{v}$ $\Theta$-Cusp Strip-Analytic Manifold (`[D, N]`, `NM-1..3`)
Defining the defect $W := (\frac{1}{2}+D)\widetilde H + \frac{1}{2}Q\widetilde G$ (compatibility $\iff W(0,s) \equiv 0$) and conjugating $W = \mathcal{N} h$ ($\frac{1}{2}Q(\mathcal{N}h) = -\mathcal{N}(uh)$, $u = w^{-1/2}$, $N(m) = \frac{\cos(\pi m)\Gamma(2m+2)}{\pi}$):
- The formal backward $O(\varepsilon^n)$ Tricomi series diverges factorially ($\sim 4^n n!$) because it is the $u\to\infty$ asymptotic expansion of the convergent forward Duhamel integral $h(u,s) = e^{-us}h_0(u) + \int_0^s e^{-u(s-s')}f(u,s')\,ds'$.
- Even with $h_0 \equiv 0$ (unmodified linear cusp $\eta_0$), $h(u,s) \sim \pi F(0,s)/u$ has its $m = \frac{1}{2}$ pole cancelled by $N(\frac{1}{2}) = 0$, enforcing $W(0,s) \equiv 0$ and $a(\tau) + D(\tau) \equiv 0$ automatically for all $s > 0$.
- In $\mathcal{X}_\rho := H^2_\rho(\frac{1}{4}) \cap H^2_\rho(\sigma_2)$ ($\sigma_2 \in (\frac{1}{2}, 1)$, e.g. $\frac{5}{8}$), the conjugated Duhamel map contracts in the Nishida–Ovsyannikov cone norm on $[0, \rho_0/a)$ with $\rho_0/a \sim c/\varepsilon > s_\infty$, yielding global weak-edge compatible solutions on $\tau \in [0,\infty)$ for all small strip-analytic $G_0 \in \mathcal{X}_{\rho_0}$.

### Claim 3 — Pure-Step ($\eta \equiv 0$) Strip-Analytic Manifold (`[D, N]`, `PM-1..3`)
Setting $\widetilde H \equiv 0$ in $(\dagger)$ gives $\widetilde G_s - \frac{1}{2}Q\widetilde G = \mathcal{P}(\widetilde G) := -(Q\widetilde G)D\widetilde G + \frac{1}{2}(Q\widetilde G)\widetilde G$, conjugated to $h(u,s) = e^{-us}h_0(u) + \int_0^s e^{-u(s-s')}(\mathcal{N}_g^{-1}\mathcal{P}(\mathcal{N}_g h))(u)\,ds'$.
- For any small strip-analytic flat seed $h_0$ (e.g. $h_0(u) = u^3 e^{-u}$), $f_{\mathcal{P}}(u,s) \sim -\frac{\pi}{2}a(\mathcal{P})/u \implies h(u,s) \sim -\frac{\pi}{2}a(\mathcal{P})/u^2 = O(u^{-2})$.
- Hence $a(\tau) \equiv 0$ ($B(\tau) \equiv A(\tau)^2 = 1$ with **zero $\Theta$-cusp**) and $c_1(\widetilde G(s)) = +3a(\mathcal{P}(\widetilde G(s))) = -\frac{3}{\pi}J(\widetilde G(s))$ hold globally and fully nonlinearly via a single-variable Nishida–Ovsyannikov contraction.

### Claim 4 — Finite-Sobolev ($H^3_{1/4}$, Non-Analytic) Pure-Step Energy Closure (`[D, N]`, `FS-1..3`)
For $A_b := \mathcal{N}_g^{-1}(b(\sigma + \partial_z))\mathcal{N}_g$ on $H^k_\sigma$ ($\sigma = \frac{1}{4}$, $m(\xi) = N(\frac{1}{4} + i\xi) = A(\xi)e^{i\Phi(\xi)}$, $W_k(\xi) = (1+\xi^2)^{k/2}$):
- In $\frac{1}{2}(A_b + A_b^*)$, the Stirling log-phase $e^{i(\Phi(\xi-\eta)-\Phi(\xi))}$ factors out identically with modulus $1$, giving $|\operatorname{Re}\langle A_b v, v\rangle_{H^k}| \le C_k\|b\|_{H^2}\|v\|_{H^k}^2$ for $k \in \{0,1,2,3\}$ (zero derivative loss).
- With hierarchical weights $c_m = 4^{3-m}$, $\sum_{m=0}^3 c_m\operatorname{Re}\langle \partial_z^m(u v), \partial_z^m v\rangle \ge \frac{1}{2}\sum_{m=0}^3 c_m\int u|\partial_z^m v|^2\,dz \ge 0$.
- Together these close a finite-Sobolev pure-step energy estimate at $k = 3$ ($h \in H^3_{1/4} \iff \widetilde G \in H^1_{1/4}$) without strip-analyticity.

---

## 3. High-Priority Attack Surfaces for Adversarial Review

Reviewers should specifically try to break the following points:

1. **Track A ($\sigma_2$ contour & $s = 0^+$ jump):**
   - In `findings_session_2026_11_21_nonlinear_cusp_manifold.md` §3 (and the self-audit note in `findings_session_2026_11_21_pure_step_nonlinear_manifold.md` §1), does the Lipschitz/splitting argument for $\mathcal{I}_0(u,s)$ in $H^2_\rho(\sigma_2)$ ($\sigma_2 = \frac{5}{8}$) actually close in the Nishida cone norm without losing a derivative on $c_0'(s)$?
   - Does the $s = 0^+$ jump in $u h(u,s)$ when $h_0 \equiv 0$ break distributional validity of $(\star)$ at $\tau = 0$, or is $h_0(u) = \frac{\pi c_0(0)}{u+1}$ required for classical/mild regularity up to $\tau = 0$?
2. **Track C (closing $\|b\|_{H^2}$ on $z \le 0$ and the rank-1 $m = \frac{1}{2}$ residue in finite Sobolev):**
   - `FS-1` bounds $|\operatorname{Re}\langle A_b v, v\rangle_{H^3}| \le C\|b\|_{H^2}\|v\|_{H^3}^2$, where $\|b\|_{H^2}$ is unweighted $H^2(\mathbb{R})$ of $b = Q\widetilde G$. Does $h \in H^3_{1/4}$ + $u$-dissipation actually control $\|b\|_{H^2(-\infty, 0]}$, or is an explicit second weight $\sigma_2 > \frac{1}{2}$ (with rank-1 residue splitting of $f_{\mathcal{P}}$ at $m = \frac{1}{2}$) strictly required? Check whether the rank-1 splitting closes cleanly in $H^3_{1/4} \cap H^4_{3/4}$ (including the $s$-regularity of $a(\mathcal{P}(\widetilde G(s)))$).
   - Also check the lower-order term $\frac{1}{2}\mathcal{N}_g^{-1}(b\widetilde G)$ in $H^3_{1/4}$: is it bounded by $C\|b\|_{H^2}\|h\|_{H^3_{1/4}}$ without phase-cancellation trouble?
3. **Distributional / weak-edge equivalence across $X = 1$:**
   - For Track B & C (pure-step $\eta \equiv 0$), verify that $a(\tau) \equiv 0$ + $G(v,\tau) = c_1(\tau)v + O(v^2)$ + $c_l = 2 - V_G$, $c_\omega = -c_l/2$ truly gives zero $\delta(X-1)$ defect AND zero $(X-1)^{-1/2}$ defect in the weak divergence form of `(eq:dyn)` on all of $\mathbb{R}$.
   - For Track A ($\sqrt{v}$ $\Theta$-cusp), verify that $\Theta = (\frac{\pi}{2} + D(\tau)\sqrt{v} + O(v^{3/2}))\mathbf{1}_{v>0}$ with $a(\tau) + D(\tau) \equiv 0$ introduces no hidden singular product or distributional mismatch at $X = 1$.
4. **Numerical gate integrity (`NM-1..3`, `PM-1..3`, `FS-1..3`):**
   - Check whether the spectral filter (`|xi| <= 0.6 xi_max` or `2/3 xi_max`) or periodic truncation on $z \in [-L/2, L/2]$ masks any real instability or edge-moment error.
   - Check `FS-3` at higher $H^k$ norms ($H^2, H^3$) on a larger domain $z \in [-48, 48]$ to confirm whether the $H^3$ drift seen in scratch runs was purely periodic-boundary wraparound from $w_g f_{\mathcal{P}} \sim e^{+z/4}$ or a real interior roughness growth.
   - Independently check $a(\tau)$ by direct real-space quadrature $-\frac{1}{\pi}\int e^{-z/2} G\,dz$ (not only via $Q(G)(z_{\mathrm{edge}})$ from the conjugated operator!).
5. **Physical-time limit $s \to s_\infty^-$ ($\tau \to +\infty$):**
   - Verify whether $V_{\widetilde G(\cdot, s)}(r^2)$ is actually Lipschitz in $r$ near $r = 0^+$ (note the $\ln(1 + r^2/w)$ factor!), so that $dr/ds = -1 + \frac{1}{2}V_{\widetilde G}(r^2)$ has a unique $C^1$ zero crossing at $s_\infty$ and $\tau(s) = \int_0^s ds'/r(s')$ diverges cleanly as $-\ln(s_\infty - s)$.
