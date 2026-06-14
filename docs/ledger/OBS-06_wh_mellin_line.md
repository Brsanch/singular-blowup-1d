---
id: OBS-06
title: Mellin-line structure of the eigenvalue problem (continuum, out-of-strip-valid) — lattice/power-series solutions are edge-divergent (exp(c/√ξ), c>0 λ-independent) for ALL λ, so eigenfunctions are Mellin-line (Ĝ_p(s)=Γ(s+λ)Γ(s+2)τ(s)); the symbol winding W=+3/2 (Re λ>0)/+½ (Re λ<0), jump +1 at Re λ=0, is the indicial-root TAUTOLOGY (not the Fredholm index: R(0)=0 on the Plancherel line Re s=0 ⟹ non-Fredholm; Coburn blocked).
category: OBS
status: verified
confidence: 0.85
dependencies: []
supports: [LEM-06, LEM-07]
cross_citations:
  - {label: "Coburn lemma", ref: "min(dim ker, dim ker*)=0 for nonzero Toeplitz symbol", note: "the index⟹no-kernel upgrade — BLOCKED here by R(0)=0 on the contour"}
verification:
  method: numerics
  rounds: 2
  adversarial_vectors:
    - "lattice-divergence reproduced (mpmath dps=40-80) — edge blow-up exp(c/√ξ), c λ-independent, both signs of Re λ ⟹ does NOT discriminate stable/unstable"
    - "winding CORRECTED: W=+3/2 (Re λ>0)/+½ (Re λ<0), jump EXACTLY +1 (per-factor exact decomposition); the earlier −0.5/jump-2 were naive-unwrap artifacts from the on-line s=0 tan-zero"
    - "R(0)=0 on Re s=0 for every λ ⟹ operator non-Fredholm on the Plancherel line; W is the raw symbol winding (3/2≠0), NOT the Fredholm index; W is the tautology Re(−λ)<0 ⟺ Re λ>0 — spectrum-independent (firewall clean) but for that reason no no-eigenvalue content"
  independent_rederivation: true
  machine_verification: "validate/edge_wh_structure.py (reproduced + corrected in workflow wf_2d8a4427)"
lean_ref: none
tags: [track2-evidence, tautology-not-certificate, mellin-line]
---

The direct run at the Wiener–Hopf factorization. It RULED OUT the lattice/termination shortcut
(eigenfunctions are Mellin-line), gave the explicit Γ-form, and CORRECTED a winding overclaim:
the winding is the indicial-root tautology, not a no-eigenvalue certificate. The genuine spectral
content is the unconstructed tan-Barnes τ(s). Filed with the index≠kernel and non-Fredholm-on-line
caveats front-and-center (asserting "W constant ⟹ no eigenvalue" would be an overclaim by exactly
the Coburn gap). Status of the leap UNCHANGED: `reduced-to-sharp-WH-factorization`. The sharpened
remaining object (construct τ(s) + side-choice + Coburn-substitute/Pincherle) is in LEM-06 and
`docs/desk/edge_mellin_symbol.md` §5.2.
