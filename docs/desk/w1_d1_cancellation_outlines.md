# W1-D1 — cancellation / continuity outlines (reopen attempt)

**Date:** 2026-08-07  
**Status:** PATH EXHAUST DONE — see `docs/findings_session_2026_08_07_w1_d1_path_exhaust.md`.  
Does not claim LEM-06 / Conjecture 2.4.  
**Parent:** `w1_plemelj_stirling_writeup.md` §3.5; magnet #4 park; LEAP_QUEUE §3.

**Wall (exact).** For \(\sigma\in(0,\tfrac12)\), \(\delta=\operatorname{Re}\lambda+2\sigma>0\),
the even-ladder functional
\[
D_1[\mathcal{F}]
=\sum_{k\ge 0}
\frac{\mathcal{F}(\tfrac12+k)}{X(\tfrac12+k)\,(\lambda+1+2k)}
\]
extends continuously to \(L^2(\operatorname{Re}s=\sigma)\) with
\(\|D_1\|=O(1/\delta)\), **without** absolute termwise / bare line density.

---

## 0. What is already ruled out (do not re-try)

| Route | Why dead |
|---|---|
| CS against bare \(k(\tau)=1/(|X||\lambda+2s|)\) on the line | Lem Def-bare⊥: \(k\sim e^{3\pi|\tau|}\) |
| \(\sum_n\sup|\Pi_n\|\) termwise | Edge depth: \(e^{2e^{|z_0|/4}}\) cancellation; LEAP_QUEUE |
| Cite Schur \(S_{\mathrm{off}}\) / KEY | Bounds dressed integral op, not scalar \(D(1)\) |
| Cite OBS-07 gains | Forcing-family forecast, not \(\|\cdot\|_{L^2\to\mathbb{C}}\) |
| Naive \(\pm\delta\) at \(s=\tfrac12\) | Analytic garbage at depth (dps-stable) |

---

## 1. Structural facts that any outline must use

1. **Real-ladder coefficients are tiny and alternating.**  
   \(a_k=1/[X(\tfrac12+k)(\lambda+1+2k)]\): \(|a_0|\sim 10^{-1}\), \(|a_1|\sim 10^{-3}\),
   then super-factorial decay; \(\operatorname{sign}(a_k)=(-1)^k\) at \(\lambda=0.3\).  
   *Termwise on the real half-integers is easy* — the hard part is defining
   \(\mathcal{F}(\tfrac12+k)\) for general line data.

2. **Plancherel line is full \(L^2\), not Hardy \(H^2\).**  
   \(L^2_\sigma\ni f \longleftrightarrow \mathcal{F}\in L^2(\sigma+i\mathbb{R})\) (bilateral).  
   Point evaluation at \(\operatorname{Re}s>\sigma\) is **not** continuous on full \(L^2\).
   Cauchy/H² formulae only see the half-line-supported (Hardy) subspace.

3. **Dense class where \(D_1\) is classical.**  
   \(\mathcal{F}_{z_0}(s)=2K_s(1)\,e^{-s z_0}\) (Laplace of \(e^{-\cosh(z-z_0)}\)) is entire;
   \(D_1(\mathcal{F}_{z_0})\) and \(\|\mathcal{F}_{z_0}\|_{L^2(\mathrm{line})}\) are well-defined.
   Probe: `validate/edge_d1_dual_probe.py`.

4. **Lem D-id.** On that dense class (mod R1), \(D(1)=C_0(\tfrac12)=D_1[\mathcal{F}]\)
   (even sum only). Bounding \(D_1\) *is* bounding the defect scalar.

---

## 2. Candidate outlines (ranked)

### Outline A — Hardy dual  (status: KILLED for \(L^2_\sigma\) bilateral data)

**Idea.** On \(H^2(\operatorname{Re}s>\sigma)\),
\[
\mathcal{F}(s_k)=\frac1{2\pi}\int_{\mathbb{R}}
\frac{\mathcal{F}(\sigma+i\tau)}{\sigma+i\tau-s_k}\,d\tau,
\quad s_k=\tfrac12+k,
\]
hence
\[
D_1[\mathcal{F}]=\int \mathcal{F}(\sigma+i\tau)\,\psi(\tau)\,d\tau,
\quad
\psi(\tau)=\frac1{2\pi}\sum_{k\ge 0}\frac{a_k}{\sigma+i\tau-s_k}.
\]
If \(\psi\in L^2(\mathbb{R})\), Riesz gives \(\|D_1\|_{H^2\to\mathbb{C}}\le\|\psi\|_2\).

**What survived.** \(\psi\in L^2(\mathbb{R})\) (probe P2: truncated norms saturate at
\(\|\psi\|_2\approx 0.122\); triangle bound matches). Early draft error: claimed
\(1/\tau\notin L^2\) at infinity — false (\(\int_{|\tau|>1}\tau^{-2}<\infty\); \(\psi\)
is regular at \(\tau=0\)). \(C_\infty=\sum a_k\neq 0\) does **not** kill \(L^2\) membership.

**What killed it for this repo.** Plancherel for \(L^2_\sigma\) is **full** \(L^2\) on
the line (bilateral \(z\in\mathbb{R}\)), not \(H^2\). On the dense class
\(\mathcal{F}=2K_s(1)e^{-s z_0}\), the H²-Cauchy formula **fails to recover**
\(\mathcal{F}(s_k)\) (rel.\ errors \(O(1)\)–\(O(2)\); findings 2026-08-07). So the
dual bound applies only on a Hardy subspace that **does not contain** the physical
dense class. Not a discharge path for D1 on \(L^2_\sigma\) without a new subspace
identification (e.g.\ two-sided Hardy split) — that would be a different outline.

**Do not revive** as “ψ∈L² ⇒ D1 bounded on L²_σ”.

### Outline B — Coboundary telescope from the line (status: **KILLED** 2026-08-07)

Targets absolute \(C_0(\tfrac12)\) via multi-line jumps. After path exhaust:
absolute \(C_0(\tfrac12)\) is **not load-bearing** for \(\mathcal{G}\) (only
differences are). Multi-line Carleson would re-solve the wrong scalar.
Detail: `findings_session_2026_08_07_w1_d1_path_exhaust.md` §3.

### Outline C — Bent Mellin–Barnes resummation (status: **KILLED** 2026-08-07)

Superseded for the resolvent by slit-plane (desk §10). Scalar Barnes for
absolute \(C_0(\tfrac12)\) is non-load-bearing after D-id retraction.
Detail: path-exhaust findings §4.

### Outline D — Slit-plane kernel for \(D(1)\) (status: **CASHED** mod Stirling 2026-08-07)

Explicit combined kernel \(K(\tau)=\kappa_{\mathrm{split}}/(X(\lambda+2s))\) with
\(|\kappa|\sim e^{-4\pi|\tau|}\) cancels \(1/|X|\sim e^{3\pi|\tau|}\) → net
\(e^{-\pi|\tau|}\) margin → \(K\in L^2\) → CS bound on slit-scalar \(D(1)\).
**Does not** restore \(D(1)=C_0(\tfrac12)\). Probe:
`validate/edge_d1_route_exhaust.py`. Detail: path-exhaust findings §2.

---

## 3. Dual probe results (2026-08-07) — decided

Script: `validate/edge_d1_dual_probe.py` + Cauchy-recovery add-on (findings).

| ID | Result |
|---|---|
| P0 | \(C_\infty\approx 0.135\neq 0\) |
| P1 | \(R\in[0.04,1.63]\) on bumps; **edge can make R large** (≠ OBS-07 annihilation) |
| P2 | \(\|\psi\|_2\) **saturates** (~0.122) — ψ∈L² as abstract H² dual |
| P3 | \(R=O(1)\) in δ on z₀=0 family — **not** clearly \(O(1/\delta)\) |
| Cauchy | H² recovery of \(\mathcal{F}(s_k)\) **fails** on bilateral bumps → A killed for L²_σ |

Detail: `docs/findings_session_2026_08_07_w1_d1_reopen.md`.

---

## 4. Discipline

- Do not green-light D1 from P1 ratios alone (family, not operator norm).
- Do not open CAP until an outline has a **written estimate** with an explicit
  constant path (magnet #4 CAP decision).
- Nonlinear / co-moving still needs its own obstruction doc before LEAP_QUEUE.

**Reopen criterion (unchanged):** a cancellation lemma *candidate with a proof
sketch of the estimate*, not another envelope table.
