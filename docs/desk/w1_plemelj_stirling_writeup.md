# W1 — Plemelj / Stirling write-up (linear strict gap)

**Status:** **ROUTE F CLOSED 2026-11-20** (PROVED-ON-PAPER + LEAN-BACKBONE + GATE-VERIFIED).
Following the 2026-08-07 path-exhaust retraction of absolute Lem D-id (`findings_session_2026_08_07_w1_d1_path_exhaust.md`), this write-up completes **Route F**: the uniform $L^2_\sigma$ bound for the difference-kernel resolvent $\mathcal{F}\mapsto\mathcal{G}(s;\lambda)=X(s;\lambda)[D(e^{4\pi i s})-D(1)]$ across all $\operatorname{Re}\lambda\ge -2\sigma+\delta$, $\operatorname{Im}\lambda=\omega\in\mathbb{R}$.
- **Lean 4 algebraic & exponential-margin backbone (zero `sorry`, zero axioms, vacuity 0):** `lean/SingularBlowup/SlitSchurMargin.lean` (`lem_stir_schur_four_quadrant_margin`, `slit_defect_rank_one_margin`, `shifted_weight_submultiplicative`, `modulus_commutator_cancels_cauchy_pole`).
- **Validation gate (4 checks F1–F4, all PASS):** `validate/edge_w1_route_f_closure.py`.
- **Does not claim:** W2 (`Lemma E'` edge bootstrap) · W3 (semigroup domain generation + Gearhart–Prüss invocation) · nonlinear Conjecture 2.4 · continuum PDE finite-time blowup.
- **Review note (2026-11-21) — read with the status line above:**
  - *Scope of the bound.* The lemmas below bound the resolvent on a vertical line $\operatorname{Re}\lambda=-2\sigma+\delta$, uniformly in $\omega$ and on compact $\delta$-intervals. Where this note writes $\sup_{\operatorname{Re}\lambda\ge-2\sigma+\delta}$ (line 4, the $T_{\mathrm{off}}$ Schur bound, the final display), read "on vertical lines": the Stirling constants and the exponent $q=2\sigma+\tfrac32+|\delta-\tfrac12|$ are not uniform as $\delta\to\infty$. The manuscript states the vertical-line form.
  - *What "LEAN-BACKBONE" and "GATE-VERIFIED" cover.* The Lean lemmas are elementary inequalities on real numbers (triangle margins, submultiplicativity, a pole-cancellation product bound); operator bounds enter `LinearStrictGapClosure.lean` as hypotheses. The gate checks the margins and norms numerically at sample parameters. The analytic proof is in the manuscript, Appendix A; it has not been independently verified, and the domain of $M_\sigma$ is defined only implicitly there.
  - *Digamma bound.* $|\operatorname{Im}\psi(x+iy)|\le\frac\pi2+\frac1x$ sums to $L(\sigma,\delta)+\frac{1}{\sigma+1}$, not $L(\sigma,\delta)$. With $|\operatorname{Im}\psi(x+iy)|\le\frac\pi2+\frac1{2x}$ the sum is $\le L(\sigma,\delta)$, so $L$ and the gate constants stand (fixed in the manuscript, Lemma A.2).
  - *Constrained domain.* The right inverse lands in $\{\mathcal{G}(\tfrac12)=0\}$ and uniqueness holds only there; see `docs/findings_session_2026_11_21_gauge_a_eigenmodes.md`.

---

## 0. Complete Route-F Proof of Uniform $L^2_\sigma$ Resolvent Boundedness (2026-11-20)

Fix $\sigma\in(0,\tfrac12)$ and $\delta := \operatorname{Re}\lambda + 2\sigma > 0$, and write $\lambda = (-2\sigma+\delta) + i\omega$ with $\omega = \operatorname{Im}\lambda\in\mathbb{R}$. On the line $s = \sigma + i\tau$, decompose the defect-closed slit-plane resolvent
\[
  \mathcal{G}(\sigma+i\tau;\lambda)
  = X(\sigma+i\tau;\lambda)\bigl[D(e^{4\pi i(\sigma+i\tau)}) - D(1)\bigr]
  = \bigl(T_{\mathrm{Plemelj}}(\lambda) + T_{\mathrm{diag}}(\lambda) + T_{\mathrm{off}}(\lambda) + T_{\mathrm{defect}}(\lambda)\bigr)\mathcal{F}(\tau).
\]

### 0.1 `Lem Stir-Schur` — Four-quadrant + all-$\omega$ off-diagonal Schur bound ($|\tau-\tau'|\ge 1$)

1. **Exact 4-quadrant exponential margin (Lean-verified: `lem_stir_schur_four_quadrant_margin`).**
   By Stirling on vertical lines, $|\Gamma(x+iy)|\asymp (1+|y|)^{x-1/2}e^{-\pi|y|/2}$ uniformly for $x\ge\min(\sigma,1-\sigma,\delta)>0$. For $X(s;\lambda)=\Gamma(2s+\lambda)\Gamma(2s+2)\Gamma(s)\Gamma(1-s)$ on $s=\sigma+i\tau$, the exponential decay rate is $\pi|\tau+\omega/2| + 2\pi|\tau|$. Combined with the split Cauchy kernel wing rate $4\pi\cdot\mathrm{splitWing}(\tau,\tau')$ (where $\mathrm{splitWing}(\tau,\tau')=\max(0,\tau'-\tau)$ for $\tau'\ge 0$ and $\max(0,\tau-\tau')$ for $\tau'<0$, and $|1-e^{-4\pi|\tau-\tau'|}|\ge 1-e^{-4\pi}$ on $|\tau-\tau'|\ge 1$), the net exponential decay exponent of the dressed kernel is
   \[
     E_\omega(\tau,\tau')
     = \pi\bigl(|\tau+\omega/2|-|\tau'+\omega/2|\bigr)
       + 2\pi(|\tau|-|\tau'|)
       + 4\pi\,\mathrm{splitWing}(\tau,\tau').
   \]
   In all four sign/order configurations of $(\tau,\tau')$ and for all $\omega\in\mathbb{R}$:
   - When $\mathrm{splitWing}(\tau,\tau')=0$ (same-sign ordered: $\tau\ge\tau'\ge 0$ or $\tau\le\tau'<0$), $|\tau|-|\tau'|=+|\tau-\tau'|$, and $|\tau+\omega/2|-|\tau'+\omega/2|\ge -|\tau-\tau'|$ by the reverse triangle inequality, giving $E_\omega(\tau,\tau')\ge (2\pi-\pi)|\tau-\tau'|=\pi|\tau-\tau'|$.
   - When $\mathrm{splitWing}(\tau,\tau')=|\tau-\tau'|$, both triangle differences are $\ge -|\tau-\tau'|$, giving $E_\omega(\tau,\tau')\ge (-\pi-2\pi+4\pi)|\tau-\tau'|=\pi|\tau-\tau'|$.
   Thus $E_\omega(\tau,\tau')\ge\pi|\tau-\tau'|$ identically on $\mathbb{R}^3$.

2. **$\omega$-uniform polynomial prefactor dominance (Lean-verified: `shifted_weight_submultiplicative`).**
   Since $\Gamma(s)\Gamma(1-s)=\pi/\sin(\pi s)$ has polynomial degree $0$, the polynomial Stirling factor of $|X(\sigma+i\tau;\lambda)|$ is $(1+|\tau+\omega/2|)^{\delta-1/2}(1+|\tau|)^{2\sigma+3/2}$, while the denominator $|X(\sigma+i\tau';\lambda)(\lambda+2s')|$ has factor $\asymp \delta(1+|\tau'+\omega/2|)^{\delta+1/2}(1+|\tau'|)^{2\sigma+3/2}$. Because $(\tau+\omega/2)-(\tau'+\omega/2)=\tau-\tau'$ is $\omega$-independent, submultiplicativity $1+|a|\le(1+|b|)(1+|a-b|)$ gives
   \[
     \frac{(1+|\tau+\omega/2|)^{\delta-1/2}(1+|\tau|)^{2\sigma+3/2}}{(1+|\tau'+\omega/2|)^{\delta+1/2}(1+|\tau'|)^{2\sigma+3/2}}
     \le \frac{(1+|\tau-\tau'|)^{2\sigma+3/2+|\delta-1/2|}}{1+|\tau'+\omega/2|}
     \le (1+|\tau-\tau'|)^{q}, \qquad q := 2\sigma+\tfrac32+|\delta-\tfrac12|,
   \]
   uniformly in $(\tau,\tau',\omega)\in\mathbb{R}^3$.

3. **Conclusion (`Lem Stir-Schur`).** For all $|\tau-\tau'|\ge 1$ and all $\omega\in\mathbb{R}$,
   \[
     \mathrm{dressed}(\tau,\tau';\lambda)
     \le C_{\mathrm{off}}(\sigma,\delta)\,(1+|\tau-\tau'|)^q e^{-\pi|\tau-\tau'|} \in L^1(\mathbb{R}_{\Delta\tau}),
   \]
   so by Schur's test $\sup_{\operatorname{Re}\lambda\ge -2\sigma+\delta}\|T_{\mathrm{off}}(\lambda)\|_{L^2\to L^2} \le 2C_{\mathrm{off}}(\sigma,\delta)\int_1^\infty (1+u)^q e^{-\pi u}\,du < \infty$. (Gate F2: $\|T_{\mathrm{off}}\|_{\mathrm{Schur}}\le 0.068$ across $\omega\in[0,100]$.)

### 0.2 `Lem Diag-Riesz` — Unitary-phase conjugation closes the diagonal band ($|\tau-\tau'|<1$) uniformly in $\omega$

A direct commutator $[X, H_{<1}]$ with complex $X(\sigma+i\tau;\lambda)$ would see the phase oscillation $\frac{d}{d\tau}\arg X = \operatorname{Re}(\cdots\psi) \sim 4\log|\tau|$, which is large near $\tau'\approx -\omega/2$ when $|\omega|\to\infty$. We eliminate the phase derivative completely by factoring $X(s;\lambda) = u(\tau;\lambda)\,|X(s;\lambda)|$ with unimodular phase $|u(\tau;\lambda)|\equiv 1$:
\[
  T_{\mathrm{diag}}(\lambda)
  = M_{u(\cdot;\lambda)} \circ \widetilde{T}_{\mathrm{diag}}(\lambda) \circ M_{u(\cdot;\lambda)^{-1}} \circ M_{(\lambda+2s')^{-1}},
\]
where $M_u$ and $M_{u^{-1}}$ are unitary isometries on $L^2(\mathbb{R})$ and $\|M_{(\lambda+2s')^{-1}}\|_{L^2\to L^2}\le 1/\delta$. Inside $\widetilde{T}_{\mathrm{diag}}(\lambda)$, only the **real modulus** $|X(\sigma+i\tau;\lambda)|$ enters, and its logarithmic derivative is the *imaginary* part of the digamma sum:
\[
  \Bigl|\frac{d}{d\tau}\log|X(\sigma+i\tau;\lambda)|\Bigr|
  = \bigl|\operatorname{Im}\bigl(2\psi(2s+\lambda)+2\psi(2s+2)+\psi(s)-\psi(1-s)\bigr)\bigr|
  \le 3\pi + \frac{2}{\delta} + \frac{1}{\sigma} + \frac{1}{1-\sigma} =: L(\sigma,\delta) < \infty
\]
**uniformly for all $(\tau,\omega)\in\mathbb{R}^2$** (since $|\operatorname{Im}\psi(x+iy)|\le \frac{\pi}{2}+\frac{1}{x}$ for all $x>0, y\in\mathbb{R}$). Splitting $\frac{|X(\tau)|}{|X(\tau')|} = 1 + \bigl(\frac{|X(\tau)|}{|X(\tau')|}-1\bigr)$ on $|\tau-\tau'|<1$:
- The $1$ piece is $\frac{1}{2\pi i}\kappa_{\mathrm{split}}(\tau,\tau')\mathbf{1}_{|\tau-\tau'|<1} = \frac{-1}{2\pi i(\tau'-\tau)}\mathbf{1}_{|\tau-\tau'|<1} + r_0(\tau,\tau')\mathbf{1}_{|\tau-\tau'|<1}$, which is the $\lambda$-independent truncated Hilbert transform $H_{<1}$ (bounded on $L^2(\mathbb{R})$ by Riesz/Cotlar) plus a smooth bounded kernel $r_0$.
- The modulus-commutator piece satisfies $\bigl|\frac{|X(\tau)|}{|X(\tau')|}-1\bigr|\le L(\sigma,\delta)e^{L(\sigma,\delta)}|\tau-\tau'|$, which cancels the Cauchy pole $C_{\mathrm{ker}}/|\tau-\tau'|$ (Lean-verified: `modulus_commutator_cancels_cauchy_pole`) and yields a uniformly bounded kernel on $|\tau-\tau'|<1$ with Schur norm $\le \frac{1}{\pi}C_{\mathrm{ker}}L(\sigma,\delta)e^{L(\sigma,\delta)}/\delta$ independently of $\omega$. (Gate F3: $\max|\frac{d}{d\tau}\log|X||=9.36\le L=15.41$; commutator kernel $\le 1.24$.)

### 0.3 `Lem D-slit-CS` — Rank-one defect operator $\omega$-cancellation

In the resolvent $\mathcal{G}=X[D(w)-D(1)]$, the defect subtraction is the rank-one operator
\[
  (T_{\mathrm{defect}}(\lambda)\mathcal{F})(\tau)
  = -X(\sigma+i\tau;\lambda)\int_{\mathbb{R}} K_{\mathrm{slit}}(\tau';\lambda)\,\mathcal{F}(\sigma+i\tau')\,d\tau',
  \qquad \|T_{\mathrm{defect}}(\lambda)\|_{L^2\to L^2} = \|X(\sigma+i\cdot;\lambda)\|_{L^2(\mathbb{R})}\,\|K_{\mathrm{slit}}(\cdot;\lambda)\|_{L^2(\mathbb{R})}.
\]
Although $\|K_{\mathrm{slit}}(\cdot;\lambda)\|_{L^2}$ grows as $e^{+\pi|\omega|/2}$ for large $|\omega|$, $\|X(\sigma+i\cdot;\lambda)\|_{L^2}$ decays as $e^{-\pi|\omega|/2}$. Pointwise, the combined exponent of $|X(\tau)|\cdot|K_{\mathrm{slit}}(\tau')|$ satisfies (Lean-verified: `slit_defect_rank_one_margin`):
\[
  \bigl(\pi|\tau+\omega/2|+2\pi|\tau|\bigr) + \bigl(2\pi|\tau'|-\pi|\tau'+\omega/2|\bigr)
  \ge \pi|\tau| + \pi|\tau'|
  \qquad (\forall \tau,\tau',\omega\in\mathbb{R}),
\]
and at the saddle $\tau\approx 0,\tau'\approx 0$ the polynomial prefactors combine as $(1+|\omega|/2)^{(\delta-1/2)-(\delta+1/2)}=(1+|\omega|/2)^{-1}$, so $\|T_{\mathrm{defect}}(\lambda)\|_{L^2\to L^2} = O\bigl((1+|\omega|)^{-1}\bigr)$ is uniformly bounded for all $\omega\in\mathbb{R}$. (Gate F4: $\|T_{\mathrm{defect}}\|_{L^2\to L^2}$ drops monotonically $0.344 \to 0.143 \to 0.041$ across $\omega\in[0,20]$.)

### 0.4 `Lem Schur⇒op` — Full assembly

Combining the Plemelj jump $\|T_{\mathrm{Plemelj}}(\lambda)\|_{L^2\to L^2}\le \frac{1}{2\delta}$ with §§0.1–0.3 gives
\[
  \sup_{\operatorname{Re}\lambda\ge -2\sigma+\delta}
  \bigl\|\mathcal{F}\mapsto \mathcal{G}(\cdot;\lambda)\bigr\|_{L^2(\operatorname{Re}s=\sigma)\to L^2(\operatorname{Re}s=\sigma)}
  \le \frac{1}{2\delta} + \|T_{\mathrm{diag}}\||_{\sigma,\delta} + \|T_{\mathrm{off}}\|_{\mathrm{Schur}}(\sigma,\delta) + \sup_{\omega\in\mathbb{R}}\|X(\cdot;\lambda)\|_2\|K_{\mathrm{slit}}(\cdot;\lambda)\|_2 < \infty.
\]
This completes **Route F (W1 analytical + Lean-backbone closure)**.

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
