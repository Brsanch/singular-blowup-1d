# W1 — Plemelj / Stirling write-up (linear strict gap)

**Status:** DRAFT — Phase 1 of magnet #4 (`research_direction_blowup_cap_2026_07_25.md`)  
**Date started:** 2026-07-25 · **Adversarial patch:** 2026-07-25 evening  
**Does not claim:** LEM-06 closed · W2/W3 · INTERNAL-CAP-CERTIFIED · continuum PDE theorem · D1 closed

**Tier labels used below:** NUMERIC-FORECAST (gate) · PROVED-ON-PAPER (target) · NOT CLAIMED / OPEN.

**Reviewer:** two-agent adversarial pass on defect bookkeeping (2026-07-25) —
skeptic FAIL on prior “default (ii) because Schur closes” non sequitur; advocate:
`D(1)≡C₀(½)` makes Route (ii) coherent. Patches below are the reconciliation.
Not a full 4-lens workflow.

---

## 0. Target statement (from LEAP_QUEUE §3 / desk §10)

For `σ ∈ (0,½)` and `δ > 0`,

\[
\sup_{\operatorname{Re}\lambda \ge -2\sigma+\delta}
\bigl\|(\lambda - M)^{-1}\bigr\|_{L^2_\sigma\to L^2_\sigma} < \infty
\]

on the defect-closed ladder realization of the convective-weight operator
(quotient by \(\mathcal{N}\)).

**Honest status after defect-growth correction (2026-07-22/25):**
representation + dressed off-diagonal Schur are a **skeleton** (de-risked).
The rank-one defect bound (Claim D1) is **OPEN** — the old decaying-kernel
story was false, so W1 is not “ink on a finished idea”; it must prove a new
estimate. Soften any “build, not a leap / no missing idea” language until D1
has a classical outline with estimates (this file §3).

---

## 1. Representation (skeleton → write-up checklist)

Under \(w = e^{4\pi i s}\), the strip \([\sigma,\sigma+\tfrac12]\) maps to
\(\mathbb{C}^*\) slit along \(\arg w = 4\pi\sigma\). Coboundary
\(C(s)-C(s+\tfrac12)=\Phi(s)\) with \(\Phi=\mathcal{F}/[X\cdot(\lambda+2\cdot)]\)
is an additive Plemelj jump across the ray.

**Claim R1 (NUMERIC-FORECAST → target PROVED-ON-PAPER).**  
The split-regularized Cauchy transform \(D\) (plain kernel on the \(w'\to 0\)
wing; one subtraction on the \(\infty\) wing) solves the jump, and Laurent
rigidity gives the defect-closed formula

\[
\mathcal{G}(s) = X(s)\bigl[D(e^{4\pi i s}) - D(1)\bigr].
\]

**Evidence already banked (relabel):** gate
`validate/edge_slitplane_representation.py` P1 — representation ≡ ladder to
rel \(10^{-22}\) on strip-interior points for Schwartz-class test
\(\widehat F\). This is an **identity on a dense/test class**, not a proof that
\(D(1)\) extends as a bounded functional on all of \(L^2_\sigma\) (that is D1).

**Write-up obligations (still open):**

1. **Lem R1 (Muskhelishvili).** Class-\(h\) hypotheses matching measured
   \(\Phi\)-growth; wing subtraction unique for restoring
   \(e^{-4\pi(\tau-\tau')}\) on the \(\infty\) wing.
2. **Growth reconciliation.** Reconcile desk \(|\Phi|\sim e^{5\pi|\tau'|/2}\)
   vs bare \(k\sim e^{3\pi|\tau|}\) (findings / `γ_def`) under the actual
   \(\tau'\)-parametrization and \(|\lambda+2s|\); state the class-\(h\)
   exponent once.
3. **Lem R2 (Laurent).** Slit-plane growth \(\Rightarrow\)
   \(\mathcal{G}=X[D(w)-D(1)]\).
4. Dense-class operator identity, then extension-by-density **after** D1.

---

## 2. Off-diagonal Schur (corrected constants)

**Claim S1 (NUMERIC-FORECAST, gate-green).**  
Dressed kernel satisfies, on the probed grid,

\[
\sup_{|\Delta\tau|\ge 1}
\bigl(\text{dressed}\bigr)\, e^{\pi|\Delta\tau|}
= S_{\mathrm{off}}(\sigma,\delta) < \infty,
\]

with effective decay \(\alpha_{\mathrm{eff}} \ge \pi\) (margin
\(\alpha_{\mathrm{eff}}-\pi \in [0.32,0.65]\) on the former corner).

**Pinned KEY line** (`validate/edge_schur_envelopes.py`, dual-run bit-identical,
anchor \(\sigma=0.4\), \(\delta=1.1\)) — **NUMERIC-FORECAST ceilings only**:

```
S_off_anchor=0.551085  alpha_anchor=3.71162  C_term_anchor=1.12545
M_diag_anchor=10.2032  defect_growth_anchor=8.89354
```

**False friend:** `defect_growth_anchor` is the log-growth rate of the *bad*
bare density \(1/(|X||\lambda+2s|)\), not a defect-control constant.
**False friend:** \(C_{\mathrm{term}}\) is a termwise-triangle proxy; LEAP_QUEUE §3
says termwise must fail at edge depth — do not cite as \(\|R\|\).

**Write-up obligations:**

1. **Lem Stir-Schur.** Prove \(\alpha_{\mathrm{eff}} \ge \pi\) from Stirling of
   \(X\) + split kernel — not from the grid alone.
2. **Lem Schur⇒op** (after D1 + Diag-Riesz): convert envelopes into
   \(\|R(\lambda)\|\lesssim C/\delta\). Until then, \(C_{\mathrm{term}}\le 1.4\)
   is shape hint only.

---

## 3. Defect functional — OPEN (load-bearing)

### 3.1 What was wrong (**Lem Def-bare⊥**)

Desk §10 previously claimed the rank-one defect kernel decays as
\(e^{-3\pi|\tau'|/2}\). Gate `edge_schur_envelopes.py` P5 shows bare density

\[
k(\tau)=\frac{1}{|X(\sigma+i\tau)|\,|\lambda+2s|}
\]

**grows** as \(e^{\gamma|\tau|}\) with \(\gamma\in[8.84,9.06]\approx 3\pi\).
Decaying-prose is **false** for \(1/|X|\) (struck from desk §10).  
**P5 PASS confirms the obstruction pin, not a bound.**

### 3.2 Claim D1 — single concrete target (OPEN)

**Claim D1 (target PROVED-ON-PAPER; currently OPEN / NOT CLAIMED).**  
For \(\sigma\in(0,\tfrac12)\), \(\delta=\operatorname{Re}\lambda+2\sigma>0\),

\[
\bigl|\mathcal{F}\mapsto D(1)[\mathcal{F}]\bigr|
= \bigl|\mathcal{F}\mapsto C_0(\tfrac12)\bigr|
\le C(\sigma,\delta)\,\|\mathcal{F}\|_{L^2(\operatorname{Re}s=\sigma)},
\]

with shape \(C\lesssim C_*/\delta\) matching bulk \(C_{\mathrm{term}}\) only as a
forecast ceiling — **not** operator-norm certified.

### 3.3 Route identity — why (ii) is the object (not a non sequitur)

**Lem D-id (to prove first).** On the dense class, with the exact even-\(n\)
evaluation at \(\tfrac12\) (desk §6: no \(\pm\varepsilon\) limit bug),

\[
D(1)=C_0(\tfrac12)
=\sum_{n\text{ even}\ge 0}
\frac{\mathcal{F}(\tfrac12+\tfrac n2)}{X(\tfrac12+\tfrac n2)\,(\lambda+1+n)}.
\]

So the Cauchy defect and the ladder sum are **the same functional**. Route (ii)
= bound that sum (or the split integral *after* using the coboundary so
cancellation precedes absolutes). Route (i) tooling (\(w\)-Jacobian / split
measure) is allowed only **inside** the cancelled integrand estimate — not as a
claim that bare \(k\in L^2\).

**Why not “default (ii) because Schur closes”:** that rationale is a **non
sequitur** (adversarial skeptic). Schur bounds the dressed integral operator;
D1 is the independent rank-one piece \(X\cdot D(1)\). Green \(S_{\mathrm{off}}\)
without D1 leaves a hole; D1 without Schur leaves the integral operator open.

**Why not bare density / triangle:** Lem Def-bare⊥ + LEAP termwise-no-go.
OBS-07 cancellation is NUMERIC-FORECAST shape evidence that *some* cancelled
object is small on test bumps — **not** a proof of D1 and **not** citable as
discharge.

### 3.4 Write-up order for D1

1. Write **Lem Def-bare⊥** (one paragraph + Stirling).
2. Write **Lem D-id** (even-\(n\) formula; named \(\pm\delta\)-limit pitfall).
3. Write **Lem Def-ladder (D1)** on dense class by cancelled ladder sum /
   cancelled split integral; extension-by-density separate.
4. Only then assemble with Schur + Diag-Riesz → op-norm sketch.

---

## 4. Diagonal band

**Claim B1-diag / Lem Diag-Riesz (classical outline).**  
On \(|\tau-\tau'|\le 1\), weight ratio \(\le e^{3\pi+}\); reduce to local
Riesz projection + \(\tfrac12\Phi\) Plemelj term. Pin against
\(M_{\mathrm{diag}}=10.203\) (NUMERIC-FORECAST) until the analytic constant
is written.

---

## 5. Circularity firewall

- Do not infer \(B=H(g)(1)=0\) from benign numerics — domain condition (Lemma D).
- Schur envelopes use explicit \(X\) and kernel only; spectrum never enters.
- **KEY / `defect_growth` ≠ defect bound** (growth of the bad density).
- **\(C_{\mathrm{term}}\) ≠ \(\|R\|\)** (triangle proxy; fails at depth).
- **Slitplane P1 ≠ arbitrary-\(L^2\) boundedness** of \(D(1)\).
- **\(S_{\mathrm{off}}\) ≠ D1** (dressed off-diagonal ≠ scalar defect).
- OBS-07/08 remain NUMERIC-FORECAST; never underwrite D1 or “skeleton complete.”
- This file is a write-up scaffold, not a proof.

---

## 6. Named lemmas (pin list)

| Name | Role |
|---|---|
| Lem R1 | Muskhelishvili class-\(h\) + wing subtraction |
| Lem R2 | Laurent ⇒ \(\mathcal{G}=X[D-D(1)]\) |
| Lem D-id | \(D(1)=C_0(\tfrac12)\) even-\(n\) identity |
| Lem Stir-X | Stirling asymptotics of \(X\) (line + rising rungs) |
| Lem Stir-Schur | Stirling + split ker ⇒ \(\alpha_{\mathrm{eff}}\ge\pi\) |
| Lem Def-bare⊥ | Bare \(k\) grows \(\gamma\approx 3\pi\); not an \(L^2\) kernel |
| Lem Def-ladder (D1) | \(\lvert\mathcal{F}\mapsto D(1)\rvert\le C(\sigma,\delta)\,\|\mathcal{F}\|_2\) |
| Lem Diag-Riesz | Diagonal band |
| Lem Schur⇒op | \(S_{\mathrm{off}}\) + Diag + D1 ⇒ \(\|R\|\lesssim C/\delta\) |

---

## 7. CAP decision (magnet #4)

| Stage | CAP needed? | Notes |
|---|---|---|
| Representation identity | Analytic first; gate = forecast on test class | Interval lift only if R1 stalls |
| Off-diagonal Schur | Analytic Stirling preferred | Route B if only grid proof survives |
| Defect (D1) | **Route unset until Lem D-id written**; then analytic Lem Def-ladder via \(D(1)\equiv C_0(\tfrac12)\) | CAP only if that estimate stalls after R1+Stir-X |
| W2/W3 | Later | — |

**Concrete choice so far:** analytic Route A; pursue D1 as Lem Def-ladder after
Lem D-id. Do **not** label INTERNAL-CAP-CERTIFIED until all four CAP items hold.

---

## 8. Session checklist (pick-up)

- [x] Desk §10 defect prose corrected (2026-07-25)
- [x] Scaffold + KEY constants
- [x] Adversarial patch: remove Schur⇒D1 non sequitur; formalize D1; Lem D-id; firewall false friends (2026-07-25)
- [ ] Growth reconciliation (\(|\Phi|\) vs \(\gamma_{\mathrm{def}}\))
- [ ] Lem R1 + Lem R2
- [ ] Lem Def-bare⊥ + Lem D-id **before** any positive D1 estimate
- [ ] Lem Def-ladder (D1) on dense class
- [ ] Lem Stir-Schur ⇒ \(\alpha_{\mathrm{eff}}\ge\pi\)
- [ ] Lem Diag-Riesz
- [ ] Do not cite bent-contour Barnes as D1 route (superseded by slit-plane §10)
- [ ] Re-run Schur + slitplane gates after any code edit
- [ ] Freeze findings when D1 cashes or re-walls as leap

**Reproduce:**

```bash
cd singular-blowup-1d
python3 validate/edge_schur_envelopes.py
python3 validate/edge_slitplane_representation.py
```
