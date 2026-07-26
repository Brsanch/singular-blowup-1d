# W1 — Plemelj / Stirling write-up (linear strict gap)

**Status:** DRAFT — Phase 1 of magnet #4 (`research_direction_blowup_cap_2026_07_25.md`)  
**Date started:** 2026-07-25  
**Does not claim:** LEM-06 closed · W2/W3 · INTERNAL-CAP-CERTIFIED · continuum PDE theorem

**Tier labels used below:** NUMERIC-FORECAST (gate) · PROVED-ON-PAPER (target) · NOT CLAIMED.

---

## 0. Target statement (from LEAP_QUEUE §3 / desk §10)

For `σ ∈ (0,½)` and `δ > 0`,

\[
\sup_{\operatorname{Re}\lambda \ge -2\sigma+\delta}
\bigl\|(\lambda - M)^{-1}\bigr\|_{L^2_\sigma\to L^2_\sigma} < \infty
\]

on the defect-closed ladder realization of the convective-weight operator
(quotient by \(\mathcal{N}\)). Skeleton: slit-plane Plemelj representation +
dressed two-weight Schur. Remaining for W1: classical write-up with
**corrected** defect bookkeeping.

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

**Evidence already banked:** gate `validate/edge_slitplane_representation.py`
— representation ≡ ladder to rel \(10^{-22}\) on strip-interior test points;
consumes only line values of \(\mathcal{F}\).

**Write-up obligations (still open):**

1. Cite Muskhelishvili class-\(h\) hypotheses matching the measured density
   growth \(|\Phi|\sim e^{5\pi|\tau'|/2}\) from \(1/X\).
2. Justify the wing-dependent subtraction as the unique choice that restores
   the \(e^{-4\pi(\tau-\tau')}\) factor on the former-corner wing.
3. State the distributional operator identity for \(D\)-defined solutions on
   dense \(L^2_\sigma\) data, then close by density.

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
anchor \(\sigma=0.4\), \(\delta=1.1\)):

```
S_off_anchor=0.551085  alpha_anchor=3.71162  C_term_anchor=1.12545
M_diag_anchor=10.2032  defect_growth_anchor=8.89354
```

**Write-up obligations:**

1. Prove \(\alpha_{\mathrm{eff}} \ge \pi\) from Stirling asymptotics of \(X\) and
   the split-kernel definition — not from the grid alone.
2. Convert \(S_{\mathrm{off}}\) + diagonal-band Riesz into an operator-norm
   bound \(\|R(\lambda)\|\lesssim C/\delta\) with explicit \(C\) (termwise proxy
   currently \(C_{\mathrm{term}}\le 1.4\) on the table — NOT YET operator norm).

---

## 3. Defect functional — CORRECTED BOOKKEEPING (load-bearing)

### 3.1 What was wrong

Desk §10 previously claimed the rank-one defect kernel
\(\mathcal{F}\mapsto C_0(\tfrac12)\) decays as \(e^{-3\pi|\tau'|/2}\).  
Gate `edge_schur_envelopes.py` shows bare density

\[
k(\tau)=\frac{1}{|X(\sigma+i\tau)|\,|\lambda+2s|}
\]

**grows** as \(e^{\gamma|\tau|}\) with \(\gamma\in[8.84,9.06]\approx 3\pi\)
(Stirling: \(|X|\sim e^{-3\pi|\tau|}\)). That decaying-prose claim is **false**
for \(1/|X|\) and has been struck from desk §10 (2026-07-25).

### 3.2 What W1 must prove instead

**Claim D1 (target PROVED-ON-PAPER).**  
The defect functional is \(L^2\)-bounded by one of:

- **(i)** a different functional kernel ( \(w\)-plane Jacobian / split measure /
  restricted \(\widehat{F}\) class ) that *is* square-integrable; or
- **(ii)** a **ladder-sum** bound that exploits alternating cancellation
  already measured at edge depth (termwise triangle must fail —
  \(e^{2e^{|z_0|/4}}\)-scale cancellation; LEAP_QUEUE §3).

**NOT CLAIMED:** that (i) or (ii) is closed in this draft. Choosing (ii) as the
default route because the slit-plane formula already encodes the ladder, and
the Schur envelopes close the dressed off-diagonal without needing bare-density
decay.

**Next write-up step:** derive an \(L^1\)/\(L^2\) bound on
\(\tau\mapsto D(1)[\mathcal{F}]\) from the subtracted Cauchy transform directly
(using only line values), with constant controlled by \(\delta=\operatorname{Re}\lambda+2\sigma\).

---

## 4. Diagonal band

**Claim B1-diag (classical outline).**  
On \(|\tau-\tau'|\le 1\), weight ratio \(\le e^{3\pi+}\); reduce to local
Riesz projection + \(\tfrac12\Phi\) Plemelj term. Standard; no new idea.
Pin the constant against \(M_{\mathrm{diag}}\) from the Schur gate
(\(M_{\mathrm{diag}}=10.203\) at anchor) as a NUMERIC-FORECAST ceiling until
the analytic Riesz constant is written.

---

## 5. Circularity firewall (inherited)

- Do not infer \(B=H(g)(1)=0\) from benign numerics — it is the generator’s
  domain condition (Lemma D).
- Schur envelopes use explicit \(X\) and kernel only; spectrum never enters.
- OBS-08 / coupled resolvent scans remain NUMERIC-FORECAST.
- This file is a write-up scaffold, not a proof.

---

## 6. CAP decision (magnet #4)

| Stage | CAP needed? | Notes |
|---|---|---|
| Representation identity | No if analytic Muskhelishvili write-up closes; float gate stays forecast | Interval lift = Route B only if write-up stalls |
| Off-diagonal Schur constants | Prefer analytic Stirling; gate constants are de-risking | If only grid proof exists → open Route B under 4-item bar |
| Defect ladder-sum | Analytic first | Same |
| W2/W3 | Later phases | — |

**Concrete choice so far:** start analytic (Route A). Do **not** label any
blowup stage INTERNAL-CAP-CERTIFIED until all four CAP items hold.

---

## 7. Session checklist (pick-up)

- [x] Desk §10 defect prose corrected (2026-07-25)
- [x] This scaffold created with pinned KEY constants
- [ ] Write §1 Muskhelishvili hypotheses + wing subtraction uniqueness
- [ ] Write §2 Stirling ⇒ \(\alpha_{\mathrm{eff}}\ge\pi\)
- [ ] Write §3.2 defect bound via ladder sum / \(D(1)\)
- [ ] Write §4 diagonal Riesz constant
- [ ] Re-run `python3 validate/edge_schur_envelopes.py` + slitplane gate after any code edit
- [ ] Freeze findings when B1 cashes or re-walls

**Reproduce (no code edit this session):**

```bash
cd singular-blowup-1d
python3 validate/edge_schur_envelopes.py
python3 validate/edge_slitplane_representation.py
```
