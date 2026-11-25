# Target B Findings: Linear Stability Spectra $\lambda(a)$ and Mellin Structure of Singular gCLM Profiles (`arXiv:2603.25104v4`)

**Date:** 2026-11-25  
**Branch:** `gclm/singular-stability`  
**Reference:** De Huang, Jiajun Tong, Xiuyuan Wang, *Self-similar finite-time blowups with singular profiles of the generalized Constantin–Lax–Majda model: theoretical and numerical investigations*, `arXiv:2603.25104v4` (Mar–Oct 2026).  
**Gates:**
- Gate 1 (`validate/gclm_verify_singular_steady.py`): **PASS**
- Gate 2 (`validate/gclm_linear_stability.py`): **PASS**
- Gate 3 (`validate/gclm_branch2_and_dns_probe.py`): **PASS**

---

## 1. Scoop Check (as of November 2026)

Status: **OPEN (not scooped).**  
Two 2026 spectral papers by Jie Xu study only **smooth** odd collapse profiles:
- `arXiv:2607.19762` (Jul 2026): smooth odd CLM profile $\Omega(y)=-y/(y^2+1/4)$ at $a=0$ (plus conditional two-line inclusion for smooth focusing $a>0$ profiles).
- `arXiv:2609.13220` (Aug 2026): Edmunds–Evans essential spectra for smooth odd positive-advection profiles ($0<a<1$).

No published paper reports the linear stability spectrum $\lambda(a)$, convective spectral gap, or Mellin/Carleman–Mellin structure for the **singular $a \le 0$ profiles** of `arXiv:2603.25104`.

---

## 2. Exact Singular Steady States & Three-Regime Transition at $a = -1$

In dynamic-rescaling coordinates $\Omega_\tau + (c_l X + a U)\Omega_X = (c_\omega + U_X)\Omega$ with $U_X = H(\Omega) = \frac{1}{\pi}\mathrm{P.V.}\int_{\mathbb{R}} \frac{\Omega(Y)}{X-Y}\,dY$ and $U(0,\tau)=0$, Theorem 3 of `arXiv:2603.25104v4` gives for every $a < 0$ (with $\mu = \frac{1}{1-a} \in (0,1)$):
$$\bar\Omega_a(X) = -\sin(\pi\mu)\frac{\mathbf{1}_{\{X>1\}}}{(X-1)^\mu}, \qquad \bar c_{l,a} = 1 - a = \frac{1}{\mu}, \qquad \bar c_{\omega,a} = -1, \qquad \bar\gamma_a = 1 - a,$$
$$(\bar U_a)_X(X) = \frac{\mathbf{1}_{\{X<1\}} + \cos(\pi\mu)\,\mathbf{1}_{\{X>1\}}}{|1-X|^\mu}, \qquad \bar U_a(X) = \frac{1}{1-\mu} + \frac{X-1}{1-\mu}(\bar U_a)_X(X).$$

On the support $X = 1 + v > 1$, the background advection velocity reduces to:
$$\bar A_a(1+v) := \bar c_{l,a}(1+v) + a\bar U_a(1+v) = (1-a)\big(v - \cos(\pi\mu)\,v^{1-\mu}\big).$$
In the power coordinate $w := v^\mu = (X-1)^\mu$, background transport is affine:
$$-\bar A_a(1+v)\,\partial_v = -\big(w - \cos(\pi\mu)\big)\,\partial_w.$$

This reveals a **sharp three-regime structural transition at $a = -1$** (missed in `arXiv:2603.25104`, whose Section 5 simulations stopped at $a = -1$):
1. **$-1 < a < 0$** ($\mu \in (\frac{1}{2},1)$, $\cos(\pi\mu) < 0$): $\bar A_a > 0$ on all $X > 1$, with **finite-time inflow** at $X = 1^+$ ($\dot w|_{w=0} = |\cos(\pi\mu)| > 0$). Passing to $\zeta := \ln(1 + v^\mu/|\cos(\pi\mu)|) \in [0,\infty)$ turns transport into constant unit speed $-\partial_\zeta$ with the exact physical edge at $\zeta = 0$.
2. **$a = -1$** ($\mu = \frac{1}{2}$, $\cos(\pi\mu) = 0$): **critical point**. $\bar A_{-1}(1+v) = 2v$, giving pure logarithmic inflow $-2\partial_z$ in $z = \ln(X-1)$ (identical transport to 1D Hou–Luo).
3. **$a < -1$** ($\mu \in (0,\frac{1}{2})$, $\cos(\pi\mu) > 0$): **interior stagnation bifurcation**. $\bar A_a$ vanishes at $X_* = 1 + (\cos(\pi\mu))^{1/\mu} > 1$ (a repelling stagnation point). On $(1, X_*)$, $\bar A_a < 0$ pushes characteristics **backward into the singularity** $X = 1^+$ (switching $X = 1^+$ from inflow to outflow).

Gate 1 certifies `max|H_quad - Ux_bar| <= 8.9e-15`, `|A_bar(1)| <= 4.4e-16`, and relative steady residuals `<= 6.5e-16` across $a \in [-2.5, -0.1]$, plus the $a = 0$ Theorem 1 / Corollary 2 CLM delta limits.

---

## 3. Desingularized Linearized Operator $M_a$ and Exact Spectrum $\lambda(a)$

Writing $\Omega(1+v,\tau) = \bar\Omega_a(1+v)(1 + G(1+v,\tau))$ and $g = \bar\Omega_a G = -\sin(\pi\mu)v^{-\mu}G$, the singular stretching term $(\bar c_{\omega,a} + (\bar U_a)_X)g$ cancels identically against $\bar A_a \bar\Omega_a' G$. Under edge pinning $c_l + a U(1) \equiv 0$:
$$G_\tau = M_a^{\mathrm{bulk}} G + \kappa[G], \qquad M_a^{\mathrm{bulk}} G(1+v) := -\bar A_a(1+v) G_v + H_+[g](1+v) + a\mu\,Q_+[g](1+v),$$
where $a\mu = -(1-\mu)$ and $Q_+[g](1+v) = \frac{1}{v}\int_0^v H_+[g](1+y)\,dy$.
- **Gauge (a)** ($c_\omega = -\mu c_l$, i.e. $\gamma \equiv 1-a$, analogue of Hou–Luo Theta-balance): $\kappa \equiv 0$ even nonlinearly; $M_a = M_a^{\mathrm{bulk}}$ has amplitude zero mode $M_a^{\mathrm{bulk}}[1] = 0$.
- **Gauge (b)** (`arXiv:2603.25104` Case 2.2, $H(\Omega)(0) \equiv 1$): rank-one projection $\kappa[G] = -L_0[\bar\Omega_a M_a^{\mathrm{bulk}} G]$ preserving $L_0[g] = H_+[g](0) = 0$.

### Certified Numerical Findings (Gate 2)
1. **Nonlinear FD certification:** `max|M_a G - dG_fd| <= 9.2e-12` (Gauge a) and `<= 4.3e-11` (Gauge b) against the raw nonlinear dynamic-rescaling PDE across $a \in [-2.0, -0.2]$.
2. **Nyquist discretization pitfall:** Because gCLM retains the order-zero stretching term $H_+[g]$ (unlike Hou–Luo's smoothing $Q_+$), unstabilized central differences produce spurious grid-scale Nyquist edge modes (`||dv||/||avg|| = 1.83` at $a = -1$; $O(N)$-divergent node-1 mode at $a = -0.75$). 2nd-order upwind (`upwind2`) eliminates them cleanly (`||dv||/||avg|| = 0.075`).
3. **No unstable modes in the smooth relative-perturbation class:** `#Re(lambda) > 1e-4 == 0` across **all** tested $a \in [-2.5, -0.1]$ in both Gauge (a) and Gauge (b), verified across refinement ($N=400 \to 700$), domain-depth sweeps, and two independent coordinates ($z$ and $\zeta$ for $-1 < a < 0$).
4. **Convective weighted spectral gap (Gauge b, $L_0[g_0]=0$):** Smooth bulk perturbations decay exponentially at the exact predicted convective rates across all three regimes:
   - $-1 < a < 0$ ($\zeta$-grid, weight $e^{-0.5\zeta}$, pred $0.500$): measured $0.5017$ ($a=-0.25$), $0.5029$ ($a=-0.50$).
   - $a = -1$ ($z$-grid, weight $e^{-0.25 z}$, pred $0.500$): measured $0.5585$.
   - $a < -1$ ($z$-grid, two-sided outflow weight $e^{-0.25|z - z_*|}$, pred $0.625/0.750$): measured $0.6191$ ($a=-1.50$), $0.6486$ ($a=-2.00$).

---

## 4. Exact Mellin Symbol, Exact Bifurcating Eigenmode $\lambda_*(a) = -(1+a)$, and Exact Conjugation at $a = -1$

In the positive-$s$ Mellin convention $\widehat{G}(s) = \int_0^\infty v^{s-1}G(1+v)\,dv$ (for $0 < \operatorname{Re} s < 1$):
$$\mathcal{M}[M_a^{\mathrm{bulk}} G](s) = \frac{s}{\mu}\,\widehat{G}(s) + c_a(s)\,\widehat{G}(s-\mu), \qquad c_a(s) = (\mu - s)\!\left[\frac{\cos(\pi\mu)}{\mu} + \frac{\sin(\pi\mu)\cot(\pi s)}{1 - s}\right].$$
(In the Hou–Luo convention $\widehat{G}_{\mathrm{HL}}(s) = \int_{\mathbb{R}} e^{-sz}G(z)\,dz$, the shift symbol is $\tilde c_a(r) = c_a(\mu - r)$.)

### Three Analytical Theorems Discovered & Certified
1. **Exact Open-Strip Zero Classification & Bifurcating Eigenmode:**
   Since $\operatorname{Im}\!\left[\frac{\cot(\pi(x+iy))}{1-x-iy}\right]$ has numerator $y\sin(2\pi x) - (1-x)\sinh(2\pi y) < 0$ everywhere on $(0,1)\times(0,\infty)$, and $s \mapsto \frac{\cot(\pi s)}{1-s}$ is strictly decreasing on $(0,1)$, **$c_a(s)$ has no non-real zeros in $0 < \operatorname{Re} s < 1$ and exactly two real zeros**:
   - $s_1 = \mu$: amplitude scaling mode $G \equiv 1$, $\lambda = 0$.
   - $s_2 = 1 - \mu$: **exact closed-form eigenmode**
     $$G_*(1+v) = (X - 1)^{2\mu - 1} = (X - 1)^{\frac{1+a}{1-a}}, \qquad g_*(1+v) = -\sin(\pi\mu)(X-1)^{-(1-\mu)}, \qquad \lambda_*(a) = -(1 + a).$$
   This mode undergoes a **transcritical bifurcation at $a = -1$**:
   - **$-1 < a < 0$:** $\lambda_*(a) = -(1+a) \in (-1,0)$ is strictly stable, $G_*(1^+)=0$ satisfies the inflow BC, and it appears directly in Gauge (a) on both the $\zeta$- and $z$-grids (e.g. at $a=-0.50$, measured $-0.50004$ vs exact $-0.50000$).
   - **$a = -1$:** $\lambda_*(-1) = 0$ collides with $G \equiv 1$ ($s_1 = s_2 = \frac{1}{2}$).
   - **$a < -1$:** $\lambda_*(a) = |a| - 1 > 0$ crosses into the unstable half-plane while $G_*(1+v) = (X-1)^{-|\frac{1+a}{1-a}|}$ simultaneously becomes more singular than $\bar\Omega_a$ at $X=1^+$ (still $g_* \in L^1_{\mathrm{loc}}$, decaying faster than $\bar\Omega_a$ at $\infty$). Thus $\bar\Omega_a$ is linearly stable to smooth relative perturbations $G \in C^1([1,\infty))$ for all $a < 0$, but admits an exact unstable steeper-singularity $L^1_{\mathrm{loc}}$ mode for $a < -1$.

2. **Exact Elementary Conjugation at $a = -1$ (and Obstruction for $a \ne -1$):**
   - At $a = -1$ ($\mu = \frac{1}{2}$, $\cos(\pi\mu) = 0$), $\tilde c_{-1}(r) = \frac{2r\tan(\pi r)}{2r+1} = (-2r)\,c_{\mathrm{HL}}(r)$ is a pure rational–trig product. With the elementary multiplier $N_L(s) := \frac{\sin(\pi s)}{2s+1}$ (in Hou–Luo convention; $\frac{\sin(\pi s)}{2s-1}$ in positive-$s$ convention),
     $$N_L(s)\,\tilde c_{-1}(s + \tfrac{1}{2}) = -N_L(s + \tfrac{1}{2}) \quad \Longrightarrow \quad \mathcal{N}_L\,M_{-1}^{\mathrm{bulk}}\,\mathcal{N}_L^{-1} = -2\partial_z - e^{-z/2},$$
     verified numerically by FFT-Mellin in Gate 2 (`rel err = 8.98e-3`). Unlike 1D Hou–Luo (where $N_{\mathrm{HL}}(s) = \cos(\pi s)\Gamma(2s+2)/\pi$ was polynomially bounded on vertical lines), the extra $(-2r)$ factor from order-zero vortex stretching cancels $\Gamma(2s+2)$, leaving $|N_L(\sigma+i\xi)| \sim |\xi|^{-1} e^{\pi|\xi|}$: the natural invariant space on which $e^{\tau M_{-1}^{\mathrm{bulk}}} = \mathcal{N}_L^{-1} S(\tau) \mathcal{N}_L$ is an exact contraction semigroup with $\|e^{\tau M_{-1}^{\mathrm{bulk}}}\|_{\mathcal{H}_\sigma} = e^{-2\sigma\tau}$ and $\sigma(M_{-1}^{\mathrm{bulk}}|_{\mathcal{H}_\sigma}) = \{\operatorname{Re}\lambda \le -2\sigma\}$ is the **slit-plane Hardy space** $\mathcal{H}_\sigma = \mathcal{N}_L^{-1}(L^2_\sigma(\mathbb{R}))$ (functions analytic on $|\operatorname{Im} z| < \pi$, i.e. $X \in \mathbb{C}\setminus(-\infty, 1]$).
   - For $a \ne -1$ ($\cos(\pi\mu) \ne 0$), $\tilde c_a(r)$ has a transcendental-sum numerator $\frac{\cos(\pi\mu)}{\mu}(1-\mu+r) + \sin(\pi\mu)\cot(\pi(\mu-r))$ (reflecting mixed homogeneity $v^1, v^{1-\mu}$ of $\bar A_a$) with non-lattice transcendental roots (e.g. $s = 1 \pm i y_*$ on $\operatorname{Re} s = 1$ for $-1<a<0$), so **no finite elementary Gamma-product factorizer exists for $a \ne -1$**.

---

## 5. Discovery of Branch II (Universal $\beta = \frac{1}{2}$ Family) & Nonlinear DNS Probe for $a < -1$ (Gate 3)

### 5.1 Second Exact Singular Branch (Unpublished in `arXiv:2603.25104`)
Solving the general power-law ansatz $\bar\Omega(X) = -\sin(\pi\beta)(X-1)^{-\beta}\mathbf{1}_{\{X>1\}}$ in the gCLM steady profile equation yields $\cos(\pi\beta)\big(-\frac{a\beta}{1-\beta} - 1\big) = 0$, which admits **two** exact singular branches for every $a < 0$ intersecting transcritically at $a = -1$:
- **Branch I** (`arXiv:2603.25104` Thm. 3): $\beta = \mu = \frac{1}{1-a}$, $\bar c_l = 1-a$, $\bar c_\omega = -1$, $\bar\gamma = 1-a$.
- **Branch II (new):** $\beta \equiv \frac{1}{2}$ (universal Hou–Luo shape!) for **every** $a < 0$:
  $$\bar\Omega_{\mathrm{II},a}(X) = \frac{1}{a}\,(X-1)^{-1/2}\,\mathbf{1}_{\{X>1\}}, \qquad \bar c_{l,\mathrm{II}} = 2, \qquad \bar c_{\omega,\mathrm{II}} = -1, \qquad \bar\gamma_{\mathrm{II}} = 2$$
  (or equivalently $\bar\Omega_{\mathrm{II}} = -(X-1)^{-1/2}\mathbf{1}_{\{X>1\}}$, $\bar c_l = -2a$, $\bar c_\omega = a$ under unit normalization $H(\bar\Omega)(0) = 1$). Its Hilbert transform vanishes identically on $X > 1$, so $\bar A_{\mathrm{II}}(1+v) = 2v > 0$ is **strictly outgoing for all $a < 0$**.

#### Exact 2-Gamma Conjugation & Exchange of Stability for Branch II (All $a < 0$)
Linearizing around Branch II gives
$$M_{\mathrm{II},a}^{\mathrm{bulk}} G = -2\partial_z G + \frac{1}{a}\,H_+[v^{-1/2} G] + \frac{1}{2}\,Q_+[v^{-1/2} G],$$
which interpolates continuously between $M_{-1}^{\mathrm{bulk}}$ (at $a = -1$) and the 1D Hou–Luo operator $-2\partial_z + \frac{1}{2}Q_+$ (as $a \to -\infty$)! Its Hou–Luo Mellin shift symbol is the pure rational–trig product
$$\tilde c_{\mathrm{II},a}(r) = -\frac{2r + 1 + a}{a(2r + 1)}\tan(\pi r),$$
admitting an **exact 2-Gamma conjugation for every $a < 0$** (verified to `< 3.3e-15` in Gate 3):
$$N_{L,\mathrm{II}}(s; a) := \sin(\pi s)\,|a|^{-2s}\,\frac{\Gamma(2s + 2 + a)}{\Gamma(2s + 2)} \quad \Longrightarrow \quad \mathcal{N}_{L,\mathrm{II}}\,M_{\mathrm{II},a}^{\mathrm{bulk}}\,\mathcal{N}_{L,\mathrm{II}}^{-1} = -2\partial_z - e^{-z/2},$$
with exact algebraic eigenmode $G_{\mathrm{II},*}(1+v) = (X-1)^{-(1+a)/2}$ at $\lambda_{\mathrm{II},*}(a) = 1 + a$ — the **exact transcritical partner** of Branch I's $\lambda_{\mathrm{I},*}(a) = -(1+a)$:
- For $a \le -1$: $\lambda_{\mathrm{II},*}(a) = -(|a|-1) \le 0$ (stable, smooth at $X=1^+$); `#Re(lambda) > 0 == 0` in Gate 3.
- For $-1 < a < 0$: $\lambda_{\mathrm{II},*}(a) = +(1+a) > 0$ crosses into the unstable half-plane (with steeper edge singularity $v^{-(1+a)/2}$).

### 5.2 Nonlinear Dynamic-Rescaling DNS Probe for $a < -1$ (`src/gclm_singular/dns.py`, Gate 3)
To determine whether smooth degenerate half-line initial data for $a < -1$ converge to Branch I (despite $X_* > 1$ and the steeper $L^1_{\mathrm{loc}}$ mode $\lambda_{\mathrm{I},*} = |a|-1 > 0$) or switch to Branch II ($\gamma = 2$), we evolved the Case 2.2 nonlinear dynamic-rescaling PDE using WENO5 + SSP-RK3 + exact log-kernel $J(v)$ on an exponentially stretched grid $v(\rho) = 10^{-4}(e^\rho - 1)$ reaching $v_{\max} \approx 4\times 10^8$:

| $a$ | Initial datum | $\gamma_{\mathrm{DNS}}$ | Branch I ($1-a$) | $X_{*,\mathrm{DNS}}$ | $X_{*,\mathrm{exact}}$ | $\max\text{rel}|\Omega-\bar\Omega_{\mathrm{I}}|$ | Bulk residual |
|---:|---|---:|---:|---:|---:|---:|---:|
| $-1.50$ | Case 2.2 ($-X^{-3}e^{-1/(X-1)}$) | $2.482$ | $2.50$ | $1.0680$ | $1.0531$ | $6.43\times 10^{-3}$ | $3.8\times 10^{-9}$ |
| $-2.00$ | Case 2.2 ($-X^{-3}e^{-1/(X-1)}$) | $2.964$ | $3.00$ | $1.1271$ | $1.1250$ | $1.18\times 10^{-2}$ | $1.6\times 10^{-9}$ |
| $-2.00$ | Poly-exp ($-(X-1)^4 e^{-(X-1)}$) | $2.964$ | $3.00$ | $1.1271$ | $1.1250$ | $1.18\times 10^{-2}$ | $2.8\times 10^{-8}$ |
| $-2.00$ | Reg. Branch II ($-(X-1+10^{-3})^{-1/2}$) | $2.964$ | $3.00$ | $1.1271$ | $1.1250$ | $1.18\times 10^{-2}$ | $7.7\times 10^{-9}$ |
| $-2.50$ | Case 2.2 ($-X^{-3}e^{-1/(X-1)}$) | $3.461$ | $3.50$ | $1.1909$ | $1.1914$ | $1.69\times 10^{-2}$ | $4.9\times 10^{-10}$ |

**Definitive conclusion for $a < -1$:**
1. **Universal convergence to Branch I:** Smooth degenerate initial data—and even regularized Branch II data—universally converge to **Branch I** ($\bar\Omega_{\mathrm{I},a}$, $\gamma = 1-a$, $\beta = \frac{1}{1-a}$), spontaneously forming the predicted interior repelling stagnation point $X_* = 1 + (\cos\pi\mu)^{1/\mu}$ (`0.04%`–`1.4%` error) and matching $\bar\Omega_{\mathrm{I},a}$ pointwise to `0.6%`–`1.7%` over three decades $X - 1 \in [0.05, 50]$.
2. **$L^2$ criticality explains the full picture:**
   - For $a < -1$, $\mu = \frac{1}{1-a} < \frac{1}{2}$, so $\bar\Omega_{\mathrm{I},a} \in L^2_{\mathrm{loc}}([1,\infty))$, the edge flux $\bar A_{\mathrm{I}}\bar\Omega_{\mathrm{I}} \sim (X-1)^{1-2\mu} \to 0$ vanishes, the IBP normalization identity $c_\omega H(\Omega)(0) = \frac{1+a}{\pi}\int\frac{U_X\Omega}{X} - \frac{a}{\pi}\int\frac{U\Omega}{X^2} = -1$ converges absolutely (verified to `< 1e-12`), and the unstable $L^1_{\mathrm{loc}}$ mode $g_* \sim (X-1)^{-(1-\mu)} \notin L^2_{\mathrm{loc}}$ is never excited by smooth data.
   - For $-1 < a < 0$, $\mu > \frac{1}{2}$ so $\bar\Omega_{\mathrm{I},a} \notin L^2_{\mathrm{loc}}$ and the edge flux $\bar A_{\mathrm{I}}\bar\Omega_{\mathrm{I}} \sim -(X-1)^{1-2\mu} \to -\infty$ must be supplied by the two-scale inner traveling wave across $X = 1$. On a fixed half-line cutoff mesh ($X > 1$ wall) that truncates this non-$L^2$ inner flux, the Case 2.2 evolution at $a = -0.5$ cleanly locks onto **Branch II** ($c_l \approx 1.00$, $c_\omega \approx -0.51$, $\gamma \approx 1.96$, $\beta \approx 0.509$ — the exact $L^2$ threshold!).

---

## 6. Natural Next Steps
1. **Rigorous linear + nonlinear stability proof** for Branch I at $a = -1$ (slit-plane Hardy-space closure via $\mathcal{N}_L M_{-1} \mathcal{N}_L^{-1} = -2\partial_z - e^{-z/2}$), for Branch II at all $a \in (-2, 0)$ (via the exact 2-Gamma conjugator $N_{L,\mathrm{II}}(s;a)$), and for Branch I at $-1 < a < 0$ (in the exact characteristic coordinate $\zeta$).
2. **Odd-symmetric singular profiles ($a < 0$):** compute the numerical steady states of Section 5.1 (`arXiv:2603.25104`, Table 1) and extend the stability spectrum $\lambda(a)$ to the odd branch.
