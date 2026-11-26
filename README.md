# singular-blowup-1d

Self-similar **finite-time blowup** of 1D nonlocal active-scalar models
(Constantin–Lax–Majda / De Gregorio / generalized CLM / 1D Hou–Luo), focusing on
**singular (unbounded) self-similar profiles**, their exact Mellin stability
spectra, two-Gamma conjugations, zero-analyticity finite-Sobolev nonlinear
closures, and machine-checked **Lean 4** + Python verification.

Author: **Bryan Sanchez** (Independent Researcher, ORCID [`0009-0000-6080-5407`](https://orcid.org/0009-0000-6080-5407))

---

## Papers & Current Status (2026-11-25)

| Target | Manuscript | Status | Lean 4 (`0 sorry`) | Verification Gates |
| :--- | :--- | :--- | :--- | :--- |
| **Target A** (1D Hou–Luo) | *Stability Spectra and Nonlinear Co-Moving Bootstrap for the Singular Self-Similar Profile of the 1D Hou–Luo Model* ([PDF](paper/stability_spectra_singular_houluo.pdf), companion to Chen–Huang–Li, `arXiv:2604.01868`) | **Submitted** to *Journal of Nonlinear Science* (25 Nov 2026) | 10 modules in [`lean/SingularBlowup/`](lean/SingularBlowup/) | `validate/` suite PASS |
| **Target B** (generalized CLM) | *Stability Spectra, Exact Two-Gamma Conjugation, and Zero-Analyticity Finite-Sobolev Nonlinear Closure for the Singular Self-Similar Profiles of the Generalized CLM Model* ([PDF](paper/stability_spectra_singular_gclm.pdf), companion to Huang–Tong–Wang, `arXiv:2603.25104v4`) | **Submitted** for peer review (25 Nov 2026) | [`lean/GclmSingularBackbone.lean`](lean/GclmSingularBackbone.lean) | Gates 1–4 (`validate/gclm_*.py`) ALL PASS |

### Cryptographic Priority Timestamps
Both manuscripts, the Lean 4 backbone, and the standalone Python verification
gates are locked under dual independent **RFC 3161** trusted timestamps
(FreeTSA.org offline-verified + DigiCert TSA backup). See
[`docs/timestamps/VERIFY.md`](docs/timestamps/VERIFY.md).

---

## Key Verified Results

### Target A — 1D Hou–Luo Singular Profile (`arXiv:2604.01868`)
- Independent certification of the closed-form singular profile
  $\bar{\Omega}(X) = (X-1)^{-1/2}\mathbf{1}_{X>1}$,
  $\bar{\Theta}(X) = \frac{\pi}{2}\mathbf{1}_{X>1}$ with $(c_l, c_\omega) = (2, -1)$
  and steady residual $\approx 2 \times 10^{-16}$.
- Edge-resolved moving-frame stability spectra in $z = \ln(Y-1)$ identifying the
  spurious fixed-grid edge artifact and resolving the convective weighted
  spectral gap $2\sigma$ modulo symmetry modes.
- Exact Mellin conjugation $\mathcal{M}_\sigma = \mathcal{N}(-2\partial_z - e^{-z/2})\mathcal{N}^{-1}$
  via $N(s) = \cos(\pi s)\Gamma(2s+2)/\pi$, yielding an explicit $C_0$-semigroup
  and vertical-line resolvent bounds on the constrained domain $\widehat{G}(\tfrac{1}{2}) = 0$.
- Co-moving nonlinear bootstrap and self-consistent dynamic-rescaling relaxation
  confirmed numerically and paired with 10 `0 sorry` Lean 4.29 modules under
  `lean/SingularBlowup/`.

### Target B — Generalized CLM Singular Profiles (`arXiv:2603.25104v4`)
- **Branch I exact Mellin symbol & spectrum:** For Huang–Tong–Wang's Branch I
  profile $\bar{\Omega}_{I,a}$ ($a < 0$), the exact Mellin symbol $c_a(s)$ has
  zero non-real zeros in $0 < \operatorname{Re} s < 1$, real roots $s = \mu$
  ($\lambda = 0$) and $s = 1 - \mu$ ($\lambda_{I,*} = -(1+a)$), and zero
  unstable eigenvalues in Gauges (a, b) across $a \in [-2.5, -0.1]$.
- **New Branch II family & transcritical exchange:** Closed-form singular family
  $\bar{\Omega}_{II,a}(X) = \frac{1}{a}(X-1)^{-1/2}\mathbf{1}_{X>1}$
  ($c_l = 2, c_\omega = -1, \gamma = 2$, valid for all $a < 0$) intersecting
  Branch I transcritically at $a = -1$ and interpolating to 1D Hou–Luo as
  $a \to -\infty$.
- **Gauge (c) exact two-Gamma conjugation (Theorem 1):** Edge-amplitude pinning
  eliminates the $\delta(X-1)$ distributional defect globally and yields an exact
  two-Gamma Mellin conjugation $\mathcal{N}_{L,II}$ of $\mathcal{M}_c$ to
  $-2\partial_z - e^{-z/2}$ on $\mathcal{H}_{\sigma,a}$ with exact spectrum
  $\{\operatorname{Re}\lambda \le -2\sigma\}$.
- **Computer-free finite-Sobolev coercivity (Theorems 2 & 4):** Rational-polynomial
  $H^2_{\sigma_-} \cap H^3_{5/8}$ coercivity for Branch II across all $a \in [-2,0)$
  and 3-line $H^1_{\sigma_-} \cap H^1_{\sigma_+}$ coercivity for Branch I across
  all $a \in [-1,0)$, machine-checked in Lean 4 (`0 sorry`, standard axioms only).
- **Unconditional zero-analyticity nonlinear closure (Theorem 3 & Appendix A):**
  Global existence, uniqueness in
  $\mathcal{X}_{\mathrm{gCLM}} \cap L^2([0,\infty); H^{5/2}_{7/8}(\mathbb{R}))$,
  exponential stability, and conservative distributional finite-time blowup via
  a Mellin-compatible weighted viscous regularizer and $O(\sqrt{\nu})$ vanishing
  viscosity.
- **$L^2$-criticality selection (DNS):** WENO5/SSP-RK3 dynamic-rescaling DNS
  demonstrates Branch I selection for $a < -1$ and Branch II locking at
  $a = -0.5$ across the $L^2$-threshold $\beta = \tfrac{1}{2}$.

---

## Reproducing the Verification

### Python Gates (Target B — gCLM)
```bash
python3 validate/gclm_verify_singular_steady.py   # Gate 1: Branch I & II profiles, residuals
python3 validate/gclm_linear_stability.py         # Gate 2: Mellin symbols, spectra, Gauges (a,b,c)
python3 validate/gclm_branch2_and_dns_probe.py    # Gate 3: Branch II transcritical & WENO5 DNS
python3 validate/gclm_rigorous_closure_gate.py    # Gate 4: Coercivity, conjugation, Branch I check
```

### Lean 4 Formalization (Lean v4.29.0 + Mathlib)
```bash
cd lean
lake build                                        # Target A: 10 modules in SingularBlowup/
lake env lean GclmSingularBackbone.lean           # Target B: gCLM coercivity & Appendix A backbone
```

### RFC 3161 Timestamp Verification
```bash
shasum -a 256 -c docs/timestamps/CHECKSUMS_TARGET_B.sha256
openssl ts -verify \
  -in docs/timestamps/CHECKSUMS_TARGET_B.freetsa.tsr \
  -data docs/timestamps/CHECKSUMS_TARGET_B.sha256 \
  -CAfile docs/timestamps/freetsa_cacert.pem \
  -untrusted docs/timestamps/freetsa_tsa.crt
```

---

## Repository Layout

```text
paper/
  stability_spectra_singular_houluo.{tex,pdf}   Target A manuscript (1D Hou–Luo)
  stability_spectra_singular_gclm.{tex,pdf}     Target B manuscript (gCLM, 21 pp)
  figures/  figures_gclm/  scripts/             Publication figures & plotting scripts
  submission_houluo/  submission_gclm/          Self-contained journal & archive bundles
lean/
  SingularBlowup/                               Target A Lean 4.29 modules (10 files, 0 sorry)
  GclmSingularBackbone.lean                     Target B Lean 4.29 backbone (0 sorry)
src/
  gclm.py  houluo.py  rescaled.py  stability*   Target A numerical engines
  gclm_singular/                                Target B profiles, operators, WENO5 DNS
validate/
  gclm_*.py                                     Target B standalone Gates 1–4
  verify_singular_steady.py  edge_*.py  ...     Target A validation gates & probes
docs/
  proof_gclm_singular_stability.md              Full analytic proof notes (Target B)
  findings_gclm_singular_stability.md           Summary of verified Target B findings
  ledger/                                       Target A proof-obligation ledger
  timestamps/                                   SHA-256 manifest + FreeTSA & DigiCert RFC 3161 tokens
```
