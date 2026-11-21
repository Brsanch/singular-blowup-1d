# W3 — Semigroup generation + Gearhart–Prüss

**Status:** **W3 WITHDRAWN 2026-11-21.** Semigroup generation and the Gearhart–Prüss decay are no longer claimed; see the Correction below. (Status on 2026-11-20 read "LINEAR STRICT-GAP CHAIN CLOSED"; kept here as history only.)  
**Date:** 2026-11-20 (updated from 2026-07-25)  
**Lean 4 backbone (0 `sorry`, 0 axioms, vacuity 0):** `lean/SingularBlowup/SlitSchurMargin.lean` (`transport_coupling_frequency_pincer`).  
**Validation gates:** `validate/edge_w1_route_f_closure.py` (F1–F4 PASS), `validate/edge_w2_w3_closure.py` (W2-A, W2-B, W3-A PASS).  
**Does not claim:** nonlinear Conjecture 2.4 · continuum PDE finite-time blowup · INTERNAL-CAP-CERTIFIED.

---

## Correction (2026-11-21) — what no longer stands

Two reviews of the manuscript on 2026-11-21 (Claude; Codex) found gaps in the W3 step, and the
manuscript now says so itself (`paper/stability_spectra_singular_houluo.tex`, Remark 1).

1. **The decay of $\delta R$ is not established.** Item 2 below writes
   $\delta R(\lambda)\mathcal{F}(s)=\frac{c(s+1/2)}{\lambda+2s}\mathcal{G}(s+\tfrac12;\lambda)$. Bounding this by
   $O((1+|\omega|)^{-1})\|\mathcal{F}\|$ needs a resolvent bound from $L^2_\sigma$ into $L^2_{\sigma+1/2}$
   (the right-hand side reads $\mathcal{G}$ on the line $\operatorname{Re}s=\sigma+\tfrac12$). W1 bounds
   $L^2_\sigma\to L^2_\sigma$ only.
2. **The inverse Laplace integral is not shown to converge.** An operator norm of size $(1+|\omega|)^{-1}$ is
   in $L^2(\mathbb{R}_\omega)$ but not in $L^1$, so the "$O(\sqrt t)$" strong-continuity step does not follow,
   and Gearhart–Prüss presupposes a $C_0$-semigroup.

Also changed since this note was written:

- The W1 resolvent bound is stated on vertical lines $\operatorname{Re}\lambda=-2\sigma+\delta$, uniformly on
  compact $\delta$-intervals, not on whole half-planes $\operatorname{Re}\lambda\ge-2\sigma+\delta$.
- W2 uniqueness holds on the constrained domain $\{\mathcal{G}(\tfrac12)=0\}$ only. Without that constraint the
  Gauge-(a) operator has an explicit eigenfunction for every $\operatorname{Re}\lambda>-2\sigma$
  (`docs/findings_session_2026_11_21_gauge_a_eigenmodes.md`, `docs/ledger/OBS-09_gauge_a_eigenmodes.md`).
- The nonlinear theorem no longer uses this step: it is a direct energy estimate
  (`docs/desk/nonlinear_comoving_leg.md` §4).

The Lean lemma `transport_coupling_frequency_pincer` is a true arithmetic inequality and is unaffected; what
is withdrawn is its use to conclude semigroup generation. §0 item 2 and the last two rows of the status
table are kept as history.

---

## 0. Complete Linear Strict-Gap Chain on $L^2_\sigma$ (2026-11-20; item 2 withdrawn 2026-11-21)

On $L^2_\sigma$ ($\sigma \in (0, \tfrac12)$) with the defect-closed domain $\mathcal{D}(M_\sigma) = R_{\mathrm{slit}}(\lambda_0) L^2_\sigma \subset \{\mathcal{G}(\tfrac12) = 0\}$:

1. **Bijectivity & Uniform Resolvent Bound in Every Half-Plane $\operatorname{Re}\lambda \ge -2\sigma + \delta$ (W1 + W2):**
   - By **W1 Route F** (`w1_plemelj_stirling_writeup.md` §0), for every $\delta > 0$, the slit-plane operator $R_{\mathrm{slit}}(\lambda) : \mathcal{F} \mapsto \mathcal{G}(\cdot;\lambda) = X(\cdot;\lambda)[D(e^{4\pi i(\cdot)}) - D(1)]$ is a holomorphic family of bounded operators on $L^2_\sigma$ satisfying
     \[
       \sup_{\operatorname{Re}\lambda \ge -2\sigma + \delta} \|R_{\mathrm{slit}}(\lambda)\|_{L^2_\sigma \to L^2_\sigma} \le C(\sigma, \delta) < \infty,
     \]
     and solves $(\lambda - M_\sigma) R_{\mathrm{slit}}(\lambda)\mathcal{F} = \mathcal{F}$ with $(R_{\mathrm{slit}}(\lambda)\mathcal{F})(\tfrac12) = 0$ identically.
   - By **W2 Left-Strip Mixed Laplace Bridge** (`w2_lemma_eprime.md` §0), $\ker(\lambda - M_\sigma) = \{0\}$ for all $\operatorname{Re}\lambda > -2\sigma$. Hence $\lambda - M_\sigma : \mathcal{D}(M_\sigma) \to L^2_\sigma$ is bijective with $(\lambda - M_\sigma)^{-1} = R_{\mathrm{slit}}(\lambda)$ for all $\operatorname{Re}\lambda > -2\sigma$.

2. **Reference-Transport Resolvent Perturbation & Semigroup Construction (W3, bypassing Lumer–Phillips):**
   - Although the energy form of $M_\sigma$ has a low-frequency positive band (`Lemma C`, so Lumer–Phillips dissipativity does not apply directly), the pure transport part $A_0 := -2\partial_z$ generates the exact translation $C_0$-group $(e^{t A_0} u)(z) = u(z - 2t)$ on $L^2_\sigma$ with $\|e^{t A_0}\|_{L^2_\sigma \to L^2_\sigma} = e^{-2\sigma t}$ and resolvent $(R_0(\lambda)\mathcal{F})(s) = \mathcal{F}(s)/(\lambda + 2s)$ (the $n = 0$ ladder rung / Plemelj jump term).
   - Subtracting $R_0(\lambda)$ from $R(\lambda) = R_{\mathrm{slit}}(\lambda)$ yields the nonlocal resolvent perturbation
     \[
       \bigl(\delta R(\lambda)\mathcal{F}\bigr)(s) := \bigl((R(\lambda) - R_0(\lambda))\mathcal{F}\bigr)(s) = \frac{c(s+1/2)}{\lambda + 2s}\,\mathcal{G}(s+\tfrac12; \lambda) = X(s;\lambda)\bigl[D_+(e^{4\pi i s}) - D(1)\bigr].
     \]
   - By the **frequency-separation pincer** (Lean: `transport_coupling_frequency_pincer`): for all $\tau, \omega \in \mathbb{R}$,
     \[
       |2\tau + \omega| + |2\tau| \ge |\omega| \implies \max\bigl(|\lambda + 2s|,\, |2s + 1|\bigr) \ge \tfrac12\sqrt{\delta^2 + \omega^2},
     \]
     so the transport resonance ($\tau \approx -\omega/2$) and the coupling support of $c(s+1/2) = O(1/(1+|\tau|))$ ($\tau \approx 0$) separate linearly in $\omega = \operatorname{Im}\lambda$. Together with our W1 high-frequency bounds $\|T_{\mathrm{off}}(\mu+i\omega)\|_{\mathrm{Schur}} = O(1/(1+|\omega|))$ (Gate F2) and $\|T_{\mathrm{defect}}(\mu+i\omega)\|_{L^2\to L^2} = O(1/(1+|\omega|))$ (Gate F4), $\delta R(\mu + i\omega)$ belongs to $L^2(\mathbb{R}_\omega; \mathcal{B}(L^2_\sigma))$ on every vertical line $\mu = -2\sigma + \delta > -2\sigma$ (and decays as $O(1/|\omega|^2)$ on test states, Gate W3-A).
   - By the **Hilbert-space Paley–Wiener theorem** on $H = L^2_\sigma$, the inverse Laplace transform
     \[
       e^{t M_\sigma} F := e^{t A_0} F + \frac{1}{2\pi}\int_{-\infty}^{\infty} e^{(\mu + i\omega)t}\,\delta R(\mu + i\omega) F\,d\omega \qquad (\mu = -2\sigma + \delta)
     \]
     defines a strongly continuous semigroup on $L^2_\sigma$ generated by $M_\sigma$ (`First Resolvent Identity` $\implies$ semigroup law; $\|\delta R(\mu+i\cdot)\|_{L^2(\mathbb{R}_\omega)} < \infty \implies \|e^{t M_\sigma}F - e^{t A_0}F\|_{L^2_\sigma} = O(\sqrt{t}) \to 0$ as $t \downarrow 0$), and Gearhart–Prüss (Engel–Nagel, *One-Parameter Semigroups*, Thm V.1.11) yields the strict linear exponential decay bound:
     \[
       \|e^{t M_\sigma}\|_{L^2_\sigma \to L^2_\sigma} \le C(\sigma, \delta)\,e^{(-2\sigma + \delta)t} \qquad (\forall t \ge 0,\; \forall \delta > 0).
     \]

---

## Status Table (2026-11-20, corrected 2026-11-21)

| Step | Tier | Status |
|---|---|---|
| Uniform resolvent $\sup_{\operatorname{Re}\lambda\ge -2\sigma+\delta}\|R(\lambda)\|<\infty$ | PROVED-ON-PAPER + LEAN-BACKBONE + GATE | **CLOSED 2026-11-20** (W1 Route F: `w1_plemelj_stirling_writeup.md` §0); 2026-11-21: stated on vertical lines and compact $\delta$-strips only, not half-planes |
| Point-spectrum uniqueness $\ker(\lambda - M_\sigma) = \{0\}$ (W2) | PROVED-ON-PAPER + LEAN-BACKBONE + GATE | **CLOSED 2026-11-20** (Left-strip mixed Laplace bridge: `w2_lemma_eprime.md` §0); 2026-11-21: on the constrained domain $\{\mathcal{G}(\tfrac12)=0\}$ only |
| Reference-transport perturbation $\delta R = R - R_0$ & $C_0$-semigroup generation | PROVED-ON-PAPER + LEAN-BACKBONE + GATE | **WITHDRAWN 2026-11-21** (see Correction; read "CLOSED 2026-11-20") |
| Gearhart–Prüss strict linear gap $2\sigma$ ($\sup = 1^-$) | PROVED-ON-PAPER (mod classical Engel–Nagel V.1.11) | **WITHDRAWN 2026-11-21** (needs W3; read "CLOSED 2026-11-20") |
