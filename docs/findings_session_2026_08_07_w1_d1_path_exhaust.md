# Findings — W1-D1 path exhaust (all routes terminal)

**Date:** 2026-08-07  
**Stop reason:** **Extremely significant finding (criterion 3(iii))** plus every live outline **B/C/D** given a terminal status.  
**Does not claim:** LEM-06 · Conjecture 2.4 · PDE finite-time blowup · INTERNAL-CAP-CERTIFIED · full Schur Stirling ink closed.

**Probe:** `validate/edge_d1_route_exhaust.py` (exit 0 on S1–S3).  
**Prior:** `docs/findings_session_2026_08_07_w1_d1_reopen.md`, magnet #4 audit, `docs/desk/w1_d1_cancellation_outlines.md`.

---

## 0. Wall statement — updated

### 0.1 Wall as previously stated (verbatim target of magnet #4 / LEAP_QUEUE §3)

For \(\sigma\in(0,\tfrac12)\), \(\delta=\operatorname{Re}\lambda+2\sigma>0\), the functional
\[
D_1[\mathcal{F}]
=\sum_{k\ge 0}
\frac{\mathcal{F}(\tfrac12+k)}{X(\tfrac12+k)\,(\lambda+1+2k)}
\]
(= “\(D(1)=C_0(\tfrac12)\)” under Lem D-id) extends continuously to
\(L^2(\operatorname{Re}s=\sigma)\) with norm \(O(1/\delta)\), without absolute termwise /
bare line density.

### 0.2 Structural correction (this session)

On the dense class \(\mathcal{F}(s)=2K_s(1)e^{-s z_0}\):

| Identity | Status |
|---|---|
| \(C_0(s)-C_0(\tfrac12)=D(w(s))-D(1)\) with \(w=e^{4\pi i s}\) | **HOLDS** (rel \({\sim}10^{-22}\)–\(10^{-25}\), probe S1b) |
| \(D(1)=C_0(\tfrac12)\) (absolute Lem D-id) | **FALSE** (rel \({\sim}0.25\)–\(0.43\), probe S1a) |
| \(C_0(s)-D(w(s))\equiv C_0(\tfrac12)-D(1)=:\Delta[\mathcal{F}]\) | **HOLDS** (s-independent nonzero functional) |
| Defect-closed \(\mathcal{G}=X[C_0-C_0(\tfrac12)]=X[D(w)-D(1)]\) | **HOLDS** (probe S2) |

**Consequence.** The resolvent map \(\mathcal{F}\mapsto\mathcal{G}\) depends only on
**differences**. Absolute equality of the two defect scalars is not used by \(\mathcal{G}\)
and is false. Bounding “the” defect as a single object \(D(1)=C_0(\tfrac12)\) was a
**misstated wall**.

**Extremely significant finding (plan criterion 3(iii)):** structural identity
collapses two-stage bookkeeping so the prior W1-D1 statement is **strictly weaker /
mis-targeted** for the resolvent. Load-bearing obligation for LEM-06’s resolvent half
is \(\|\mathcal{F}\mapsto\mathcal{G}\|_{L^2\to L^2}<\infty\) via the **difference kernel**
(slit Cauchy \(D(w)-D(1)\), already the Schur skeleton), not absolute continuity of
the even ladder alone.

---

## 1. Terminal status table

| Route | Terminal | One-line reason |
|---|---|---|
| **A** H²-Cauchy dual (bilateral \(L^2_\sigma\)) | **KILLED** (prior 2026-08-07) | Cauchy recovery fails on bilateral dense class; dual only on Hardy subspace |
| Bare \(k\sim e^{3\pi|\tau|}\) / Lem Def-bare⊥ | **KILLED** (prior) | Absolute density grows; CS on bare \(1/|X|\) impossible |
| Termwise \(\sum\sup|\Pi_n|\) | **KILLED** (prior) | Edge depth alternating cancellation forbids triangle |
| Schur \(S_{\mathrm{off}}\) / KEY as D1 discharge | **KILLED** (prior) | Bounds dressed op, not absolute defect scalar; still not a D1 proof |
| OBS-07 gains as \(\|\cdot\|_{L^2\to\mathbb{C}}\) | **KILLED** (prior) | Forcing family, not operator norm |
| \(\pm\delta\) eval at \(s=\tfrac12\) | **KILLED** (prior) | Analytic garbage at depth |
| **B** Coboundary telescope | **KILLED** | Recovers absolute \(C_0(\tfrac12)\) from jumps — non-load-bearing after §0.2; multi-line Carleson would only bound the wrong scalar |
| **C** Bent Mellin–Barnes scalar | **KILLED** | Superseded for \(\mathcal{G}\) by slit-plane (desk §10); remaining Barnes work targets absolute ladder / full resolvent resummation already replaced — no independent D1 discharge path |
| **D** Slit-plane kernel for \(D(1)\) | **CASHED** (mod Stirling on combined kernel) | Explicit \(K\in L^2\) forecast: \(|\kappa_{\mathrm{split}}|\sim e^{-4\pi|\tau|}\) cancels \(1/|X|\sim e^{3\pi|\tau|}\) leaving \(e^{-\pi|\tau|}\) margin; CS \(\|D(1)\|\le\|K\|_2\|\mathcal{F}\|_2\). **Does not** restore absolute Lem D-id |
| **E** (new) Absolute even-ladder \(C_0(\tfrac12)\) continuity on full \(L^2\) | **LEAP-PARKED** | Still open as a pure ladder problem; **not load-bearing** for \(\mathcal{F}\mapsto\mathcal{G}\) once slit differences are used |
| **F** (new) Difference-kernel Schur Stirling ink (true resolvent remainder) | **LEAP-PARKED** | Exact remaining: uniform bound of \(\mathcal{F}\mapsto X[D(w)-D(1)]\) with fully written Stirling for off-diagonal \(e^{-\pi|\Delta\tau|}\) margin + diagonal Riesz band; generic gate re-runs do not write the estimate |

---

## 2. Outline D — cash sketch (hypotheses + constant path)

**Definition (line data).** With \(\Phi=\mathcal{F}/[X(\lambda+2\cdot)]\) and split factors as in
`validate/edge_slitplane_representation.py`,
\[
D(1)[\mathcal{F}]
=\frac{-1}{2\pi i}\int_{\mathbb{R}}
\kappa_{\mathrm{split}}(\tau;w{=}1)\,
\frac{\mathcal{F}(\sigma+i\tau)}{X(\sigma+i\tau)\,(\lambda+2(\sigma+i\tau))}
\,d\tau
=\int_{\mathbb{R}} K(\tau)\,\mathcal{F}(\sigma+i\tau)\,d\tau.
\]

**Elementary factors.** \(w'=e^{4\pi i\sigma}e^{-4\pi\tau}\). At \(w=1\):
- \(\tau\to+\infty\): \(|\kappa|\asymp e^{-4\pi\tau}\);
- \(\tau\to-\infty\) (subtraction wing): \(|\kappa|\asymp e^{-4\pi|\tau|}\).

**Stirling (Lem Stir-X, already cashed).** \(|X(\sigma+i\tau)|^{-1}\sim e^{3\pi|\tau|}\) times poly.

**Net.** \(|K(\tau)|\le C(\sigma,\delta)\,e^{-\pi|\tau|}(1+|\tau|)^{N}\) for explicit \(N\) from
four \(\Gamma\) prefactors + \(1/|\lambda+2s|\). Hence \(K\in L^2(\mathbb{R})\) and
\[
\bigl|D(1)[\mathcal{F}]\bigr|\le\|K\|_2\,\|\mathcal{F}\|_{L^2(\mathrm{line})}.
\]
\(\delta\)-dependence sits in \(1/|\lambda+2s|\) and \(\Gamma(2s+\lambda)\) distance to poles
(forecast \(O(1)\)–mild in \(\delta\) on bumps; do not overclaim \(O(1/\delta)\) for this scalar).

**Hypotheses.** Split Cauchy definition of \(D\) on the dense class / line data as in the
gate; Stirling as in Lem Def-bare⊥ / Stir-X. Muskhelishvili uniqueness of full \(D(w)\) as
a function (R1) is **not** required for the scalar CS bound on this explicit integral.

**NUMERIC-FORECAST:** probe S3 — \(\|K\|_2\) saturates (\({\approx}0.099\) at \(\lambda=0.3\));
\(|K|e^{\pi|\tau|}\) decays in \(|\tau|\).

**What D-cash does *not* do:** prove absolute \(C_0(\tfrac12)\) ladder continuity; prove
LEM-06; identify \(D(1)\) with \(C_0(\tfrac12)\).

---

## 3. Outline B — kill detail

Coboundary \(C_0(s)-C_0(s+\tfrac12)=\Phi(s)\) telescopes to express absolute
\(C_0(\tfrac12)\) as a sum of line values of \(\Phi\) on shifted lines after contour
motion. That targets **absolute** \(C_0(\tfrac12)\). After §0.2:

1. Absolute \(C_0(\tfrac12)\) is not the resolvent’s defect object.
2. Contour shift still requires multi-line bounds equivalent to a new Carleson/Schur
   problem for a non-load-bearing scalar.
3. No reduction of \(\mathcal{F}\mapsto\mathcal{G}\) beyond the slit difference kernel
   already available.

**Terminal: KILLED** as a D1/resolvent discharge route.

---

## 4. Outline C — kill detail

Bent Mellin–Barnes was the pre-slit plan to resum the ladder and Schur-test a Gamma
kernel (vertical rotation fails at slope \(\pi/4\)). Desk §10 **superseded** that for
the resolvent by slit-plane Plemelj + dressed Schur. Remaining “Barnes for scalar
\(C_0(\tfrac12)\)” is again absolute ladder, non-load-bearing after §0.2, and does not
revive a path that avoids either (D) or (F).

**Terminal: KILLED** as an independent D1 discharge route.

---

## 5. What remains (honest LEM-06 remainder — not D1 theater)

| Item | Status |
|---|---|
| Representation \(\mathcal{G}=X[D(w)-D(1)]\) on dense class | Machine-checked (existing P1 gate) |
| Off-diagonal dressed envelopes NUMERIC-FORECAST | Existing KEY; Lem Stir-Schur ink still open |
| Diagonal band Riesz outline | Open ink |
| Absolute Lem D-id | **Retracted** |
| Absolute ladder continuity (route E) | LEAP-PARKED, non-load-bearing |
| Lemma E′ / Gearhart–Prüss (W2/W3) | Unchanged, blocked on full resolvent bound |
| PDE blowup theorem | Out of scope |

**Exact inequality left for the true gap (route F, LEAP-PARKED):**  
For \(\sigma\in(0,\tfrac12)\), \(\delta>0\),
\[
\sup_{\operatorname{Re}\lambda\ge -2\sigma+\delta}
\bigl\|(\lambda-M)^{-1}\bigr\|_{L^2_\sigma\to L^2_\sigma}<\infty,
\]
equivalently uniform \(L^2\)-bound of \(\mathcal{F}\mapsto X[D(w(\cdot))-D(1)]\), by a
written Stirling estimate on the **difference** kernel (not absolute \(D(1)\), not
OBS-07, not bare \(1/|X|\)).

Why generic continuations fail: re-running Schur tables / dual probes only forecasts
constants; absolute-ladder attacks solve the wrong scalar; H² dual does not see
bilateral \(L^2_\sigma\).

---

## 6. False-friend integrity

This findings doc does **not** treat as a D1 or LEM-06 proof:

- Schur envelopes / \(S_{\mathrm{off}}\) / KEY  
- OBS-07 depth gains  
- bare \(k\sim e^{3\pi|\tau|}\)  
- H² dual on bilateral data  
- bump-family Rayleigh quotients alone as operator norm  

---

## 7. What was measured

| Check | Result |
|---|---|
| Abs \(D(1)\) vs \(C_0(\tfrac12)\) | rel \({\sim}0.25\)–\(0.43\) on bumps |
| Difference identity | rel \({\sim}10^{-22}\) |
| \(\|K\|_2\) saturation | yes (\({\sim}0.099\) at anchor) |
| \(|K|e^{\pi|\tau|}\) | decaying / bounded |

---

## 8. Next human leap (if any)

Write Lem Stir-Schur + diagonal band for the **difference** kernel (route F) —
classical Stirling/Riesz ink on an already-identified operator. Absolute ladder
continuity is optional curiosity, not the gap wall.

---

## 9. Artifact index

- This file (final path-exhaust)  
- `validate/edge_d1_route_exhaust.py`  
- `docs/desk/w1_d1_cancellation_outlines.md` (prior outlines)  
- `docs/findings_session_2026_08_07_w1_d1_reopen.md` (Outline A kill)  
