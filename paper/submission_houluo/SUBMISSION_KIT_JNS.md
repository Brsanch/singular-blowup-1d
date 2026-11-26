# Submission Kit — *Journal of Nonlinear Science* (Springer) & Zenodo Preprint
**Manuscript:** `paper/stability_spectra_singular_houluo.pdf`
**Author:** Bryan Sanchez (Independent Researcher)
**ORCID:** `0009-0000-6080-5407` (`https://orcid.org/0009-0000-6080-5407`)
**Email:** `bryansanchez.bas@gmail.com`

---

## 1. Portal Links
- **Zenodo Preprint Upload:** https://zenodo.org/uploads/new
- ***Journal of Nonlinear Science* Homepage & Portal:** https://www.springer.com/journal/332 (click *Submit manuscript*)

---

## 2. Copy-Paste Metadata (Editorial Manager / SNAPP & Zenodo)

### Title
```text
Stability Spectra and Nonlinear Co-Moving Bootstrap for the Singular Self-Similar Profile of the 1D Hou-Luo Model
```

### Plain-Text Abstract (234 words)
```text
We investigate the linear and nonlinear stability of the singular unbounded self-similar blowup profile of the one-dimensional Hou-Luo model, omega_t + u omega_x = theta_x, theta_t + u theta_x = 0, u_x = H[omega], constructed by Chen-Huang-Li (arXiv:2604.01868, 2026): Omega_bar(X) = (X-1)^{-1/2} 1_{X>1}, Theta_bar(X) = (pi/2) 1_{X>1}, with scaling exponents (c_l, c_omega) = (2, -1).

First, in the edge-resolving coordinate z = ln(Y-1) and convective weighted space L^2_gamma(R, e^{-2 gamma z} dz), we compute the linear stability spectra across two independent discretizations and certify #{Re lambda > 0} = 0 for support-interior relative perturbations.

Second, we establish an exact Mellin-Carleman conjugation of the linearized nonlocal operator M to the local transport-damping operator -2 d_z - e^{-z/2} via an explicit Gamma-function Fourier multiplier, characterizing the continuous spectrum {Re lambda <= -2 gamma} and proving exact semigroup contraction bounds on Hardy-weighted spaces.

Third, we establish a nonlinear co-moving bootstrap framework closing small support-interior perturbations of the singular profile in analytic/Hardy-weighted spaces, together with an explicit accounting of the edge jump and distributional structure at X = 1.

All algebraic identities and symbol bounds are machine-checked in Lean 4, and all spectral, operator, and nonlinear numerical checks are reproduced by standalone Python verification gates.
```

### MSC 2020 Codes
```text
35B44, 35B35, 35C06, 35Q35, 76B03
```

### Keywords
```text
1D Hou-Luo model; singular self-similar blowup; linear stability spectrum; Mellin transform; Carleman-Mellin conjugation; nonlinear stability; Lean 4 formalization
```

### Declarations (Required by Springer)
- **Data/Code Availability:** `Companion Python verification scripts and Lean 4 formalizations are archived in the repository https://github.com/Brsanch/singular-blowup-1d and in the accompanying Zenodo deposit.`
- **Funding:** `No funding was received for conducting this study.`
- **Competing Interests:** `The author has no competing interests to declare that are relevant to the content of this article.`

---

## 3. Cover Letter to the Editor (*Journal of Nonlinear Science*)

```text
Dear Editors of the Journal of Nonlinear Science,

Please find enclosed the manuscript "Stability Spectra and Nonlinear Co-Moving Bootstrap for the Singular Self-Similar Profile of the 1D Hou-Luo Model" for consideration as an Original Article in the Journal of Nonlinear Science.

Following the recent discovery of an explicit singular unbounded self-similar blowup profile for the 1D Hou-Luo model by Chen-Huang-Li (arXiv:2604.01868, 2026), we determine its linear stability spectrum, derive an exact Mellin-Carleman conjugation of the linearized nonlocal operator to -2 d_z - e^{-z/2}, and develop a nonlinear co-moving stability bootstrap with Lean 4 formal verification of key algebraic symbol bounds and standalone Python reproducibility gates.

This manuscript is original, has not been published previously, and is not under consideration at any other journal. (A distinct companion study analyzing the generalized Constantin-Lax-Majda model is being submitted separately.)

Sincerely,
Bryan Sanchez
Independent Researcher
ORCID: 0009-0000-6080-5407
Email: bryansanchez.bas@gmail.com
```

---

## 4. Suggested Referees
1. **De Huang** — Peking University (`dhuang@math.pku.edu.cn`)
2. **Jiajie Chen** — Courant Institute, New York University
3. **Thomas Y. Hou** — California Institute of Technology
4. **In-Jee Jeong** — Seoul National University
5. **Pavel M. Lushnikov** — University of New Mexico
