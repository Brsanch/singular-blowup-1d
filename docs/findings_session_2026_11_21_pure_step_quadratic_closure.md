# Full Time-Dependent $O(\varepsilon^2)$ Pure-Step Weak-Edge Closure (2026-11-21)

**Context:** Completes Codex's 2026-11-21 turning-point quadratic profile construction (`/private/tmp/singular-quadratic-5qaifq1t/prototype.py`) for the pure $\Theta$-step ($\eta \equiv 0$, $A \equiv B \equiv \xi \equiv 1$, $c_l = 2 - \varepsilon V_{G_1} - \varepsilon^2 V_{G_2}$, $c_\omega = -c_l/2$) 1D Hou–Luo weak-edge compatibility problem.

**Evidence tiers:** `[D]` hand derivation; `[N]` numerical checks (`validate/edge_pure_step_quadratic_closure.py`, `PQ-1..PQ-3 PASS`). Not Lean-formalized; higher orders $O(\varepsilon^3)$ and full nonlinear series convergence remain open.

---

## 1. Complete $O(\varepsilon^2)$ Solution (`[D, N]`)

Expanding $G(v,\tau) = \varepsilon G_1(v,\tau) + \varepsilon^2 G_2(v,\tau) + O(\varepsilon^3)$ around the single-seed flat linear solution $G_1(v,\tau) = (e^\tau/\beta(\tau))^3 \Phi_1(v/\beta(\tau)^2)$ ($\beta(\tau) = (b+1)e^\tau - 1$), the exact $O(\varepsilon^2)$ PDE is
$$\partial_\tau G_2 + 2v\partial_v G_2 - \tfrac{1}{2}Q G_2 = \underbrace{-(QG_1)v\partial_v G_1 + \tfrac{1}{2}(QG_1)G_1}_{S_{\mathrm{nl}}(v,\tau)} + \underbrace{V_{G_1}(\tau)\,v\partial_v G_1}_{S_V(v,\tau)}.$$

1. **Scale-covariant piece $G_{2,\mathrm{hom}}$ (Codex's turning-point profile):**
   Setting $G_{2,\mathrm{hom}}(v,\tau) := (e^\tau/\beta(\tau))^6 \psi(v/\beta(\tau)^2)$ reduces $(\partial_\tau + 2v\partial_v - \frac{1}{2}Q)G_{2,\mathrm{hom}} = S_{\mathrm{nl}}$ to $-6\psi - 2\psi_z - \frac{1}{2}Q\psi = \Phi_{1,b}(4\Phi_1 + \Phi_{1,b})$. Under conjugation $\psi = \mathcal{N}h$, this becomes the forced first-order ODE
   $$u h'(u) + (u - 6)h(u) = f(u), \qquad f = \mathcal{N}^{-1}[\Phi_{1,b}(4\Phi_1 + \Phi_{1,b})], \qquad h(6) = 0,$$
   integrated stably outward on both sides from the homogeneous turning point $u = 6$. Its leading edge coefficient is $c_1(\psi) = -\frac{3}{\pi}D_0 = -179.5284722663 = -\frac{3}{\pi}J$, cancelling $a' = -\frac{1}{3}c_1 - \frac{1}{\pi}J$ identically, while $a(\psi) = 0$.

2. **Gauge-modulation piece $G_{2,V}$ (exact closed form via conjugation):**
   In the conjugated variable $k = \mathcal{N}^{-1}G$, $(\partial_\tau - u\partial_u + u)k_{2,V} = -\frac{1}{2}V_{G_1}(\tau)e^{3\tau}(3u^3 - \beta u^4)e^{-\beta u}$ has the exact decaying solution
   $$G_{2,V}(v,\tau) = A(\tau)\,G_1(v,\tau) - \tfrac{1}{2}B(\tau)\,(Q G_1)(v,\tau),$$
   $$A(\tau) = \frac{3}{2}\int_\tau^\infty V_{G_1}(s)\,ds, \qquad B(\tau) = -\frac{1}{2}e^\tau\int_\tau^\infty e^{-s}\beta(s)V_{G_1}(s)\,ds.$$
   Both $G_1$ and $Q G_1 = 2e^{3\tau}\partial_\beta\Phi_\beta$ have identically zero $m=0$ moment ($a(G_1) \equiv 0$, $a(QG_1) = 2\partial_\beta[\widehat{\Phi}_\beta(\frac{1}{2})] \equiv 0$) and are flat at $v=0$, so $G_{2,V}$ preserves $a(G_2)(\tau) \equiv 0$ and $c_1(G_2) = -\frac{3}{\pi}J(G_1)$ for all $\tau \ge 0$.

---

## 2. Verification Summary (`validate/edge_pure_step_quadratic_closure.py`)

- **`PQ-1`:** Bulk homological residual `< 5e-8` (`n=8192`) vs `source_sup = 9.17`; $c_1(\psi) = -\frac{3}{\pi}D_0$ exact to $10^{-12}$.
- **`PQ-2`:** Restoring the exact $c_1 e^z$ edge tail $-\frac{2c_1}{\pi}e^{z_{\min}/2}$ yields $a(\psi) = -3.18\times 10^{-7}$ on $[-16,16]$ (`-1.62e-8` on $[-20,20]$ at `n=16384`).
- **`PQ-3`:** Full time-dependent $O(\varepsilon^2)$ PDE residual for $G_2 = G_{2,\mathrm{hom}} + G_{2,V}$ across $\tau \in \{0, 0.3, 0.8\}$, $v \in \{0.05, 0.3, 1.5\}$ achieves `max rel < 3e-6` (`1.0e-9` at $\tau=0, v=0.05$).

---

## 3. Scope & Open Questions

- **Established (`[D, N]`):** Even without a $\sqrt{v}$ $\Theta$-cusp, the pure-step $O(\varepsilon^2)$ quadratic obstruction $J > 0$ is cancelled globally in time by a curved correction $G_2 = G_{2,\mathrm{hom}} + G_{2,V}$ with $a(G_2)(\tau) \equiv 0$ and $c_1(G_2)(\tau) = -\frac{3}{\pi}J(G_1)(\tau)$.
- **Open:** Extending to $O(\varepsilon^3)$ and proving convergence of a full nonlinear pure-step compatible manifold (noting that $c_1(G_2) \neq 0$ introduces a nonzero $O(\varepsilon^2)$ edge slope whose $O(\varepsilon^3)$ feedback into $F_1$ and $c_2$ must be checked).
