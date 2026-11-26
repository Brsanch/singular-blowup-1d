# Submission Kit — *Nonlinearity* (IOP / LMS) & Zenodo Preprint
**Manuscript:** `paper/stability_spectra_singular_gclm.pdf`
**Author:** Bryan Sanchez (Independent Researcher)
**ORCID:** `0009-0000-6080-5407` (`https://orcid.org/0009-0000-6080-5407`)
**Email:** `bryansanchez.bas@gmail.com`

---

## 1. Portal Links
- **Zenodo Preprint Upload:** https://zenodo.org/uploads/new
- ***Nonlinearity* ScholarOne Portal:** https://mc04.manuscriptcentral.com/non
  (via https://iopscience.iop.org/journal/0951-7715 → *Submit an article*)

---

## 2. Copy-Paste Metadata (ScholarOne & Zenodo)

### Title
```text
Stability Spectra, Exact Two-Gamma Conjugation, and Zero-Analyticity Finite-Sobolev Nonlinear Closure for the Singular Self-Similar Profiles of the Generalized CLM Model
```

### Plain-Text Abstract (258 words — safely under ScholarOne's 300-word cap)
```text
We analyze the linear stability spectra, Mellin-Carleman operator structure, nonlinear selection, and rigorous linear and nonlinear stability of the singular self-similar blowup profiles of the one-dimensional generalized Constantin-Lax-Majda (gCLM) model, omega_t + a u omega_x = u_x omega, u_x = H(omega), for a < 0, extending the half-line outer singular family (Branch I) constructed by Huang-Tong-Wang (arXiv:2603.25104v4, 2026).

First, we discover a second explicit half-line singular family (Branch II), Omega_{II,a}(X) = (1/a)(X-1)^{-1/2} 1_{X>1}, valid for every a < 0 with (c_l, c_omega, gamma) = (2, -1, 2). Branch II intersects Branch I transcritically at a = -1, exchanging exact partner eigenmodes lambda_{I,*}(a) = -(1+a) and lambda_{II,*}(a) = 1+a, and interpolates to the singular 1D Hou-Luo profile as a -> -infinity.

Second, using singularity-adapted logarithmic and characteristic coordinates with second-order upwind discretization, we compute the linear stability spectra across a in [-2.5, -0.1] and certify #{Re lambda > 0} = 0 for smooth relative perturbations.

Third, introducing an edge-amplitude-pinning gauge (Gauge (c)) that eliminates the delta(X-1) distributional defect identically, we prove four theorems: (1) an exact two-Gamma Mellin conjugation to -2 d_z - e^{-z/2} with spectrum {Re lambda <= -2 sigma}; (2) computer-free rational-polynomial coercivity in H^2_{sigma_-} cap H^3_{5/8} for all a in [-2, 0); (3) unconditional global existence, uniqueness, and exponential stability in pure finite Sobolev space (zero analyticity, zero codimension-1 constraint) via weighted viscous regularization, covering Branch I at a = -1 and Branch II for a in [-2, 0); and (4) elementary H^1_{sigma_-} cap H^1_{sigma_+} coercivity of Branch I for all a in [-1, 0).

Fourth, WENO5/SSP-RK3 dynamic-rescaling DNS reveals an L^2-criticality selection mechanism across a = -1. All algebraic positivity certificates are machine-checked in Lean 4 (0 sorry), and all spectral, PDE, and DNS results are verified by four standalone Python gates.
```

### MSC 2020 Codes
```text
35B44, 35B35, 35C06, 35Q35, 76B03
```

### Keywords (7 — ScholarOne allows up to 10)
```text
generalized Constantin-Lax-Majda model; singular self-similar blowup; linear stability spectrum; Mellin transform; Carleman-Mellin conjugation; finite-Sobolev nonlinear stability; Lean 4 formalization
```

### Data & Code Availability Statement
```text
All Python verification gates (validate/gclm_*.py), source modules (src/gclm_singular/), figure scripts (paper/scripts/make_gclm_figures.py), the Lean 4.29 / Mathlib formalization (lean/GclmSingularBackbone.lean), and SHA-256 / RFC 3161 reproducibility records are openly archived in the repository https://github.com/Brsanch/singular-blowup-1d (branch gclm/singular-stability) and in the accompanying Zenodo deposit.
```

### Competing Interests & Funding Statement
```text
The author declares no competing financial or non-financial interests. This research received no external funding.
```

---

## 3. Cover Letter to the Editor (*Nonlinearity*)

```text
Dear Editors of Nonlinearity,

Please find enclosed the manuscript "Stability Spectra, Exact Two-Gamma Conjugation, and Zero-Analyticity Finite-Sobolev Nonlinear Closure for the Singular Self-Similar Profiles of the Generalized CLM Model" for consideration as a research paper in Nonlinearity.

Following the classic introduction of the generalized Constantin-Lax-Majda (gCLM) model in Nonlinearity (Okamoto-Sakajo-Wunsch, 2008) and the recent construction of explicit singular half-line self-similar blowup profiles for a < 0 by Huang-Tong-Wang (arXiv:2603.25104v4, 2026), the linear stability spectra, Mellin-Carleman operator structure, and nonlinear stability of these singular profiles remained open. In this manuscript we:

1. Uncover a second explicit half-line singular family (Branch II) valid for every a < 0, which intersects Huang-Tong-Wang's Branch I transcritically at a = -1 with an exact exchange of partner eigenmodes and interpolates continuously to the singular 1D Hou-Luo profile as a -> -infinity.
2. Compute the linear stability spectra lambda(a) across a in [-2.5, -0.1] in singularity-adapted logarithmic and characteristic coordinates, certifying #{Re lambda > 0} = 0 for smooth relative perturbations and diagnosing spurious central-difference Nyquist edge artifacts.
3. Introduce an edge-amplitude-pinning dynamic-rescaling gauge (Gauge (c)) that eliminates the delta(X-1) distributional defect identically and yields four rigorous theorems: exact two-Gamma Mellin conjugation on slit-plane Hardy space (Theorem 1), computer-free rational-polynomial H^2 cap H^3 coercivity for all a in [-2, 0) (Theorem 2), unconditional global existence, uniqueness, and exponential stability in pure finite Sobolev space requiring zero analyticity (Theorem 3 and Appendix A), and an elementary three-line H^1 coercivity proof for Huang-Tong-Wang's Branch I across all a in [-1, 0) (Theorem 4).
4. Explain nonlinear branch selection across a = -1 via an L^2-criticality edge-flux mechanism confirmed by WENO5/SSP-RK3 dynamic-rescaling DNS.

All algebraic positivity certificates and bootstrap inequalities are machine-checked in Lean 4 (zero sorry, standard axioms only), and all continuous spectral, PDE, and DNS claims are reproduced by four standalone Python verification gates.

This manuscript is original, has not been published previously, and is not under consideration for publication elsewhere. A preprint version is archived on Zenodo / GitHub in accordance with IOP Publishing's preprint policy.

Sincerely,
Bryan Sanchez
Independent Researcher
ORCID: 0009-0000-6080-5407
Email: bryansanchez.bas@gmail.com
```

---

## 4. Suggested Referees (ScholarOne asks for 3–5)
1. **De Huang** — School of Mathematical Sciences, Peking University (`dhuang@math.pku.edu.cn`) — *Singular self-similar blowup of gCLM and 1D Hou–Luo*
2. **Jiajun Tong** — Beijing International Center for Mathematical Research, Peking University — *Singular self-similar gCLM profiles (`arXiv:2603.25104`)*
3. **Tarek M. Elgindi** — Department of Mathematics, Duke University — *gCLM self-similar blowup and stability*
4. **In-Jee Jeong** — Department of Mathematical Sciences, Seoul National University — *gCLM advection vs. stretching and self-similar blowup*
5. **Pavel M. Lushnikov** — Department of Mathematics and Statistics, University of New Mexico — *gCLM collapse vs. blowup and spectral stability*

---

## 5. Exact Files to Upload
- **To Zenodo (first):**
  1. `paper/stability_spectra_singular_gclm.pdf`
  2. `paper/submission_gclm/gclm_code_and_lean_v1.1.zip`
- **To *Nonlinearity* ScholarOne (second):**
  - Upload **only** `paper/stability_spectra_singular_gclm.pdf` as *Main Document* (IOP accepts single-PDF initial submissions; choose **Subscription** for \$0 publication charge).
