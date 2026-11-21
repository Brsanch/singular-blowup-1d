# W2 — Lemma E′ (eliminated via Left-Strip Mixed Laplace Bridge)

**Status:** **CLOSED 2026-11-20** (PROVED-ON-PAPER + LEAN-BACKBONE + GATE-VERIFIED; edge bootstrap bypassed via left-strip mixed Laplace bridge)  
**Date:** 2026-11-20 (updated from 2026-07-25)  
**Lean 4 backbone (0 `sorry`, 0 axioms, vacuity 0):** `lean/SingularBlowup/SlitSchurMargin.lean` (`left_strip_edge_exponent_le`, `left_strip_far_exponent_le`, `shift_half_weight_exponent_eq`, `left_strip_gamma_real_floors`).  
**Validation gate (W2-A, W2-B, W2-C PASS):** `validate/edge_w2_w3_closure.py`.  
**Does not claim:** nonlinear Conjecture 2.4 · continuum PDE finite-time blowup · INTERNAL-CAP-CERTIFIED.  
**Review note (2026-11-21) — read with the status line above:**

- *Uniqueness holds on the constrained domain only.* The step "$KG=e^{-z/2}T_cG\in L^2_{\sigma-1/2}$" identifies the nonlocal term with $\tfrac12 Q_{\mathrm{op}}G$. For the Gauge-(a) operator $-2\partial_z+\tfrac12\widetilde{Q}_{\mathrm{op}}$ the two differ by the constant $\tfrac12H[g](1)$, whose Laplace transform on $[0,\infty)$ has a pole at $s=0$; then $P(0)\ne0$ is allowed and $P\equiv c_0\ne0$ gives a genuine eigenfunction $\mathcal{G}=c_0X$ for every $\operatorname{Re}\lambda>-2\sigma$. So the proposition is a statement about the operator restricted to $\{\mathcal{G}(\tfrac12)=0\}$ (equivalently $H[g](1)=0$, the frozen-exponent operator). Details and a numerical probe: `docs/findings_session_2026_11_21_gauge_a_eigenmodes.md`, `docs/ledger/OBS-09_gauge_a_eigenmodes.md`.
- *What "LEAN-BACKBONE" and "GATE-VERIFIED" cover.* The Lean lemmas are exponent inequalities and an integer budget; the gate checks two identities numerically for one test function. The Morera gluing across the strip boundary and the full-strip Parseval Laurent-orthogonality bound $\frac{1}{2}\sum_{k\in\mathbb{Z}}|c_k|^2 e^{-8\pi k\tau} = \int_{\sigma-1/2}^\sigma |P(\sigma'+i\tau)|^2\,d\sigma'$ (which eliminates single-ray cancellation across $k$) are written out in the manuscript (Lemma A.4) and spot-checked in `W2-C`.

---

## 0. Complete Proof of Proposition E Unconditional via the Left-Strip Mixed Laplace Bridge (2026-11-20)

**Why the right-strip edge bootstrap was a detour.** Prior formulations of Proposition E (`convective_weight_spectrum.md` §5) worked on the *right* strip $\sigma \le \operatorname{Re} s \le \sigma + \tfrac12$, which required proving a priori that an $L^2_\sigma$ eigenfunction $G$ carries extra $e^{(\sigma+1/2)z}$ edge decay as $z \to -\infty$ (`Lemma E'`). Because $(2\partial_z + \lambda)\mathcal{G}(s) = c(s+\tfrac12)\mathcal{G}(s+\tfrac12)$ relates lines separated by $\tfrac12$, **the left strip $\sigma - \tfrac12 \le \operatorname{Re} s \le \sigma$ wraps onto the exact same $w = e^{4\pi i s}$ slit plane $\mathbb{C}^* \setminus \{\arg w = 4\pi\sigma\}$**, and on the left strip both half-line embeddings hold *automatically* by elementary weight monotonicity:

1. **Edge half-line $z \in (-\infty, 0]$ on the left strip (Lean: `left_strip_edge_exponent_le`).**
   For $z \le 0$ and $\sigma' \in [\sigma - \tfrac12, \sigma]$, $-2\sigma' z \le -2\sigma z$, so $e^{-2\sigma' z} \le e^{-2\sigma z}$ and $G|_{(-\infty, 0]} \in L^2_{\sigma'}((-\infty, 0])$ for all $\sigma' \in [\sigma - \tfrac12, \sigma]$. By Paley–Wiener,
   \[
     \mathcal{G}_-(s) := \int_{-\infty}^0 G(z) e^{-sz}\,dz
   \]
   is holomorphic on $\operatorname{Re} s < \sigma$ with $\sup_{\sigma-1/2 \le \sigma' \le \sigma} \|\mathcal{G}_-(\sigma'+i\cdot)\|_{L^2(\mathbb{R})} \le \|G\|_{L^2_\sigma} < \infty$.

2. **Far-field half-line $z \in [0, \infty)$ on the left strip (Lean: `shift_half_weight_exponent_eq`, `left_strip_far_exponent_le`).**
   Because $T_c : L^2_\sigma(\mathbb{R}) \to L^2_\sigma(\mathbb{R})$ is a bounded Fourier multiplier ($\sigma \in (0, \tfrac12)$ avoids all poles of $\tan(\pi s)$) and multiplication by $e^{-z/2}$ is an exact isometry from $L^2_\sigma(\mathbb{R})$ to $L^2_{\sigma-1/2}(\mathbb{R})$ ($-z - 2(\sigma-\tfrac12)z = -2\sigma z$), every $G \in L^2_\sigma(\mathbb{R})$ has $KG = e^{-z/2}T_c G \in L^2_{\sigma-1/2}(\mathbb{R})$. On the far field $z \ge 0$, $-2\sigma' z \le -2(\sigma-\tfrac12)z$ for all $\sigma' \in [\sigma-\tfrac12, \sigma]$, so $(KG)|_{[0, \infty)} \in L^2_{\sigma'}([0, \infty))$ across the entire left strip, and
   \[
     \mathcal{F}_+(s) := \int_0^\infty (KG)(z) e^{-sz}\,dz
   \]
   is holomorphic on $\operatorname{Re} s > \sigma - \tfrac12$ with $\sup_{\sigma-1/2 \le \sigma' \le \sigma} \|\mathcal{F}_+(\sigma'+i\cdot)\|_{L^2(\mathbb{R})} \le \|KG\|_{L^2_{\sigma-1/2}} < \infty$.

3. **The mixed half-line bridge $\Phi_{\mathrm{mixed}}(s)$ (Gate W2-A: error $< 3\times 10^{-31}$).**
   For any $G \in L^2_\sigma(\mathbb{R})$ solving $(2\partial_z + \lambda)G = KG$ with $\operatorname{Re}\lambda + 2\sigma = \delta > 0$, define on the left strip $\sigma - \tfrac12 \le \operatorname{Re} s \le \sigma$:
   \[
     \Phi_{\mathrm{mixed}}(s) := (\lambda + 2s)\,\mathcal{G}_-(s) + 2G(0) + \mathcal{F}_+(s).
   \]
   By Steps 1–2, $\Phi_{\mathrm{mixed}}(s)$ is **holomorphic on the open left strip $\sigma - \tfrac12 < \operatorname{Re} s < \sigma$**.
   - On the **right boundary $\operatorname{Re} s = \sigma$**, where $G|_{[0,\infty)} \in L^2_\sigma([0,\infty))$, integration by parts on $[0,\infty)$ gives $2G(0) + \mathcal{F}_+(s) = (\lambda + 2s)\mathcal{G}_+(s)$, so $\Phi_{\mathrm{mixed}}(s) = (\lambda + 2s)\mathcal{G}(s)$.
   - On the **left boundary $\operatorname{Re} s = \sigma - \tfrac12$**, where $G|_{(-\infty,0]} \in L^2_{\sigma-1/2}((-\infty,0])$, integration by parts on $(-\infty,0]$ gives $(\lambda + 2s)\mathcal{G}_-(s) + 2G(0) = \int_{-\infty}^0 (KG)(z)e^{-sz}\,dz$, so $\Phi_{\mathrm{mixed}}(s) = \mathcal{L}[KG](s) = c(s+\tfrac12)\mathcal{G}(s+\tfrac12)$.

4. **Holomorphic $\tfrac12$-periodic quotient and automatic $s = 0$ vanishing (Lean: `left_strip_gamma_real_floors`, Gate W2-B).**
   In $(\lambda + 2s)X(s) = \Gamma(2s + \lambda + 1)\Gamma(2s + 2)\Gamma(s)\Gamma(1 - s)$, across the entire closed left strip $\sigma' \in [\sigma - \tfrac12, \sigma]$ we have $2\sigma' + \operatorname{Re}\lambda + 1 \ge \delta > 0$, $2\sigma' + 2 \ge 2\sigma + 1 > 1$, and $1 - \sigma' \ge 1 - \sigma > \tfrac12$, while $1/\Gamma(s)$ is entire with a simple zero at $s = 0 \in (\sigma - \tfrac12, \sigma)$. Hence
   \[
     P(s) := \frac{\Phi_{\mathrm{mixed}}(s)}{(\lambda + 2s)X(s)} = \frac{\Phi_{\mathrm{mixed}}(s)}{\Gamma(2s+\lambda+1)\Gamma(2s+2)\Gamma(s)\Gamma(1-s)}
   \]
   is **holomorphic on the open left strip $\sigma - \tfrac12 < \operatorname{Re} s < \sigma$**, satisfies $P(0) = 0$ (since $\Phi_{\mathrm{mixed}}(0)$ is finite and $1/\Gamma(0) = 0$), and has matching boundary traces on $\operatorname{Re} s = \sigma - \tfrac12$ and $\operatorname{Re} s = \sigma$:
   \[
     P(\sigma - \tfrac12 + i\tau)
     = \frac{c(\sigma+i\tau)\mathcal{G}(\sigma+i\tau)}{c(\sigma+i\tau)X(\sigma+i\tau)}
     = \frac{\mathcal{G}(\sigma+i\tau)}{X(\sigma+i\tau)}
     = P(\sigma + i\tau).
   \]
   Thus $P(s)$ extends to an entire $\tfrac12$-periodic function on $\mathbb{C}$ with $\mathcal{G}(\sigma+i\tau) = X(\sigma+i\tau)P(\sigma+i\tau) \in L^2(\mathbb{R}_\tau)$. By `LEM-08` (`PointSpectrum.lean`: Laurent modes must satisfy $|k| < 3/4 < 1$), $P(s) \equiv c_0$ is constant, and $P(0) = 0$ forces $c_0 = 0 \implies \mathcal{G} \equiv 0 \implies G \equiv 0$. $\blacksquare$

---
## NUMERIC-FORECAST evidence (not a proof)

Gate `validate/edge_eprime_probe.py` · findings
`docs/findings_session_2026_07_22_eprime_residual_floor.md`.

At \(\sigma=0.4\), \(\lambda=0.3\):
```
KEY gamma_w=1.10035 r_min=0.585763 lam_spur_re=1.7066
    alpha_free=-0.544389 alpha_bc=3.02852 maxRe_bc=-0.00156221
```

| Fact | Reading |
|---|---|
| \(r_{\min}=0.586\); \(r(-14)/r(-6)=55.6\) | No constructible approx eigenfunction from \(G_{\mathrm{hom}}\); residual **grows** with depth |
| \(\gamma_w=\operatorname{Re}\lambda+2\sigma\) | Homogeneous branch not in \(L^2_\sigma\) (measured) |
| Free-edge spur \(\operatorname{Re}\lambda\approx +1.71\) killed by inflow BC | Discrete BC matches the continuum object E′ names |

**Falsification probe:** a continuum counterexample to E′ would need a sequence
with residual \(\to 0\) as edge depth \(\to\infty\). Observed discrete residual
moves the wrong way — supports keeping E′, does not prove it.

---

## Why not CAP this session

E′ is a boundary-regularity / Liouville-type statement. No fixed-point
formulation with computable defect bounds is ready. Opening CAP without
conditioning would violate “condition before certifying.”

---

## Park decision

- Keep E′ as a **named hypothesis** in Prop E.
- Do not weaken Prop E to hide the gap.
- Do not cite residual-floor KEY as PROVED-ON-PAPER.
- Reopen W2 only with a continuum Mellin/Liouville outline or a CAP-ready
  formulation — not more discrete residual scans.
