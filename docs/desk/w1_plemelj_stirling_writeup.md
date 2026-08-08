# W1 — Plemelj / Stirling write-up (linear strict gap)

**Status:** PATH EXHAUST 2026-08-07 — absolute Lem D-id **retracted**; Claim D1 as
prior wall **mis-targeted**. See `findings_session_2026_08_07_w1_d1_path_exhaust.md`.
Lem Def-bare⊥ still cashed; slit \(D(1)\) CS sketch cashed; true remainder = difference-kernel Schur ink.
**Does not claim:** LEM-06 · Conjecture 2.4 · W2/W3 · INTERNAL-CAP-CERTIFIED · continuum PDE

**Tier labels:** NUMERIC-FORECAST · PROVED-ON-PAPER (mod named hyp) · OPEN / PARKED.

**Reviewer:** two-agent adversarial pass on defect bookkeeping (2026-07-25), then
solo lemma write + park. Not a 4-lens workflow.

---

## 0. Target and cashout posture

Uniform resolvent on \(L^2_\sigma\):
\(\sup_{\operatorname{Re}\lambda\ge -2\sigma+\delta}\|(\lambda-M)^{-1}\|<\infty\).

| Piece | Tier | Status |
|---|---|---|
| Lem Def-bare⊥ | PROVED-ON-PAPER | **Cashed** (§3.1) |
| Lem D-id | PROVED-ON-PAPER mod R1 dense-class identity | **Cashed** (§3.3) |
| Lem Stir-X (line growth) | PROVED-ON-PAPER (Stirling) | **Cashed** (§3.1) |
| Claim D1 / Lem Def-ladder | RETRACTED as wall | Absolute D-id false (§0 path exhaust) |
| Lem Stir-Schur / R1 full / Diag / Schur⇒op | OPEN | Dependent on D1 or independent ink |

---

## 1. Representation (unchanged obligations)

Slit-plane: \(\mathcal{G}(s)=X(s)[D(e^{4\pi i s})-D(1)]\).  
Gate P1 (`edge_slitplane_representation.py`): ≡ ladder to \(10^{-22}\) on
Schwartz test class — **NUMERIC-FORECAST** dense-class identity, not \(L^2\)
boundedness of \(D(1)\).

**Still open for full R1:** Muskhelishvili class-\(h\) write-up; growth
reconciliation \(|\Phi|\sim e^{5\pi|\tau'|/2}\) vs bare \(k\sim e^{3\pi|\tau|}\);
wing-subtraction uniqueness.

---

## 2. Off-diagonal Schur (forecast ceilings)

KEY (`edge_schur_envelopes.py`, \(\sigma=0.4\), \(\delta=1.1\)):
```
S_off=0.551085  α_eff=3.71162  C_term=1.12545  M_diag=10.2032  γ_def=8.89354
```
Lem Stir-Schur (\(\alpha_{\mathrm{eff}}\ge\pi\) analytic) remains OPEN ink.
\(C_{\mathrm{term}}\) ≠ \(\|R\|\); \(\gamma_{\mathrm{def}}\) ≠ defect bound.

---

## 3. Defect functional

### 3.1 Lem Def-bare⊥ — **PROVED-ON-PAPER**

**Statement.** Let \(s=\sigma+i\tau\), \(\sigma\in(0,\tfrac12)\),
\(X(s)=\Gamma(2s+\lambda)\Gamma(2s+2)\Gamma(s)\Gamma(1-s)\), and
\[
k(\tau)=\frac{1}{|X(\sigma+i\tau)|\,|\lambda+2s|}.
\]
Then as \(|\tau|\to\infty\), \(\log k(\tau)/|\tau|\to\gamma\) with
\(\gamma=3\pi\) at leading Stirling order (up to \(o(1)\) and polynomial-in-
\(|\tau|\) prefactors from \(\Gamma\)). In particular \(k\notin L^2(\mathbb{R})\)
and \(k\) does **not** decay as \(e^{-3\pi|\tau|/2}\).

**Proof sketch.** Standard Stirling:
\(|\Gamma(x+iy)|\sim\sqrt{2\pi}\,|y|^{x-1/2}e^{-\pi|y|/2}\) as
\(|y|\to\infty\) at fixed \(x\). The four \(\Gamma\)-factors in \(X\) on the line
\(\operatorname{Re}s=\sigma\) contribute three copies of the
\(e^{-\pi|\tau|/2}\) exponential from the arguments with imaginary part
\(\sim\tau\) or \(2\tau\) in the combination recorded in desk/LEM-09
(\(|X|\sim e^{-3\pi|\tau|}\) leading), hence \(1/|X|\) grows as
\(e^{3\pi|\tau|}\). The factor \(1/|\lambda+2s|\) is \(O(1/|\tau|)\). ∎

**NUMERIC-FORECAST confirmation:** gate P5,
\(\gamma_{\mathrm{def}}\in[8.84,9.06]\subset(2.5\pi,3.5\pi)\),
`defect_growth_anchor=8.89354`.

**Corollary.** No \(L^2\) defect bound of the form
\(|C_0(\tfrac12)|\le\int k|\mathcal{F}|\) (or CS against bare \(k\)) can hold
for general line data. Absolute termwise ladder bounds at edge depth are
likewise forbidden (OBS-07 / LEAP_QUEUE: \(e^{2e^{|z_0|/4}}\) cancellation).

### 3.2 Claim D1 — target (still OPEN)

\[
\bigl|\mathcal{F}\mapsto D(1)[\mathcal{F}]\bigr|
=\bigl|\mathcal{F}\mapsto C_0(\tfrac12)\bigr|
\le C(\sigma,\delta)\,\|\mathcal{F}\|_{L^2(\operatorname{Re}s=\sigma)}.
\]

### 3.3 Lem D-id — **PROVED-ON-PAPER modulo R1 dense-class identity**

**Statement.** On the dense class \(\mathcal{D}\) of line data for which
(i) the Prop F ladder
\[
C_0(s)=\sum_{n\ge 0}\frac{\mathcal{F}(s+n/2)}{X(s+n/2)\,(\lambda+2s+n)}
\]
converges in the strip (even-\(n\) closed form at \(s=\tfrac12\); **never**
naive eval or \(\pm\delta\) symmetric limit — desk §6), and
(ii) the slit-plane formula
\(\mathcal{G}=X[D(w)-D(1)]\) agrees with the defect-closed ladder
\(\mathcal{G}=X[C_0-C_0(\tfrac12)]\) pointwise on a set of strip-interior
points with accumulation (gate P1: rel \(10^{-22}\) on Schwartz tests),
one has
\[
D(1)=C_0(\tfrac12)
=\sum_{n\text{ even}\ge 0}
\frac{\mathcal{F}(\tfrac12+n/2)}{X(\tfrac12+n/2)\,(\lambda+1+n)}.
\]

**Proof.** Under (ii),
\[
C_0(s)-C_0(\tfrac12)=D(e^{4\pi i s})-D(1)
\quad\text{for strip-interior }s.
\]
Hence \(C_0(s)-D(e^{4\pi i s})\) is independent of \(s\). Both representations
are normalized by the same defect-closure condition \(\mathcal{G}(\tfrac12)=0\)
(Lemma D / Prop F). The unique constant implementing that normalization in the
ladder is \(C_0(\tfrac12)\); in the slit-plane formula it is \(D(1)\). Agreement
of \(\mathcal{G}\) on an accumulating set forces the constants to coincide. The
even-\(n\) sum is the only safe evaluation of \(C_0(\tfrac12)\) (desk §6). ∎

**Hypotheses not discharged here:** full Muskhelishvili uniqueness of \(D\)
(Lem R1); extension of (ii) from test class to all of \(\mathcal{D}\) by a
classical identity theorem (currently NUMERIC-FORECAST via P1).  
**Pitfall:** \(\pm\delta\) limits at \(s=\tfrac12\) are analytically wrong at
depth (dps-stable spurious gains \(10^{12}\)–\(10^{15}\)).

### 3.4 Why Route (ii) after Lem D-id

Once Lem D-id holds, bounding \(D(1)\) **is** bounding the even ladder sum.
Schur envelopes do not discharge this (adversarial 2026-07-25). Bare \(k\)
and triangle are ruled out by Lem Def-bare⊥.

### 3.5 Lem Def-ladder (D1) — **PARKED (named wall W1-D1)**

**Attempted routes and why they do not cash:**

| Route | Failure |
|---|---|
| CS / bare density | Contradicts Lem Def-bare⊥ |
| Termwise \(\sum_n\sup|\Pi_n|\) | Forbidden at depth (LEAP_QUEUE §3); \(C_{\mathrm{term}}\) is bulk proxy only |
| Cite OBS-07 gains | NUMERIC-FORECAST on a forcing family; not \(\|\cdot\|_{L^2\to\mathbb{C}}\) |
| Cite Schur \(S_{\mathrm{off}}\) | Bounds dressed integral operator, not scalar \(D(1)\) |
| Cancelled ladder / Abel | **No classical estimate written** — this is the missing idea |

**Wall statement (exact).** Prove, for \(\sigma\in(0,\tfrac12)\) and
\(\delta=\operatorname{Re}\lambda+2\sigma>0\), that the even-\(n\) functional
\[
\mathcal{F}\mapsto
\sum_{n\text{ even}\ge 0}
\frac{\mathcal{F}(\tfrac12+n/2)}{X(\tfrac12+n/2)\,(\lambda+1+n)}
\]
is continuous on \(L^2(\operatorname{Re}s=\sigma)\), with norm \(O(1/\delta)\),
by a cancellation argument that does **not** pass through
\(\sum|a_n|\,\|\mathcal{F}\|\) or bare \(k\in L^2\).

**Classification:** not technical ink on a finished skeleton — **leap-shaped
estimate** (or multi-week analytic project). Magnet #4 parks here per DONE WHEN
(“honest park if obstruction is not technical execution”). Re-queued in
LEAP_QUEUE §3.

---

## 4. Diagonal band — outline only

Lem Diag-Riesz: \(|\Delta\tau|\le 1\), weight \(\le e^{3\pi+}\), local Riesz +
\(\tfrac12\Phi\). Pin \(M_{\mathrm{diag}}=10.203\) as forecast ceiling.
**Not cashed** (blocked behind op-norm assembly needing D1).

---

## 5. Circularity firewall

Unchanged: no spectrum-from-gains; KEY/γ_def/C_term/S_off/P1 false friends;
OBS-07/08 never discharge D1.

---

## 6. Lemma status table

| Name | Tier | Status |
|---|---|---|
| Lem Def-bare⊥ | PROVED-ON-PAPER | Cashed |
| Lem D-id | PROVED-ON-PAPER mod R1 | Cashed |
| Lem Stir-X (line) | PROVED-ON-PAPER | Cashed (inside Def-bare⊥) |
| Lem Def-ladder (D1) | RETRACTED as wall | Absolute D-id false; ladder non-load-bearing (2026-08-07) |
| Lem D-slit-CS | CASHED mod Stirling | Combined kernel \(K\in L^2\); path-exhaust §2 |
| Lem Stir-Schur / Diag / difference-kernel op | OPEN / LEAP-PARK | True resolvent remainder (route F) |
| Lem R1, R2 | OPEN | Deferred |

---

## 7. CAP decision (cashed)

**No stage labeled INTERNAL-CAP-CERTIFIED.** D1 stalled before any constant
needed interval arithmetic. Concrete choice: **analytic-only for W1; do not
open HQSW-style CAP on blowup until W1-D1 has an estimate outline worth
certifying.** B4 vacuously held (no fake certificate).

Harness-generality test (§4 magnet): **not transferred** this arc — HQSW remains
the sole INTERNAL-CAP-CERTIFIED instance.

---

## 8. Checklist

- [x] Def-bare⊥ + D-id
- [x] D1 parked as W1-D1 with exact wall statement
- [x] CAP decision: analytic-only / no fake cert
- [x] 2026-08-07 reopen attempt: outlines + dual probe
      (`docs/desk/w1_d1_cancellation_outlines.md`,
      `docs/findings_session_2026_08_07_w1_d1_reopen.md`,
      `validate/edge_d1_dual_probe.py`). H² dual killed for bilateral L²_σ;
      next chip candidates: slit kernel for D(1), bent Barnes scalar, or ship paper.
- [ ] Future: reopen only with a cancellation lemma candidate (not more OBS-07)

```bash
python3 validate/edge_schur_envelopes.py
python3 validate/edge_slitplane_representation.py
```
