---
id: OBS-05
title: Argument-principle eigenvalue count of M_GG — #eig(M_GG,edge_bc=True) with Re λ≥{0.05,0.25,0.5,1} is 0 across N∈{1400,2000}, z_min∈{−12,−14,−16}, and under a g-mass penalty; the contour winding (1/2πi)∮tr((M−λ)^{-1})dλ over {0.1≤Re≤6,|Im|≤8} = 0 (Im ~1e-15). Robust IN-STRIP corroboration that no unstable Ω-shape point mode exists — but NOT the Wiener–Hopf index.
category: OBS
status: verified
confidence: 0.85
dependencies: []
supports: [LEM-06, LEM-07]
cross_citations: []
verification:
  method: numerics
  rounds: 2
  adversarial_vectors:
    - "DIFFERENT winding: (C) winds the spectral parameter λ (counts matrix eigenvalues); the WH index winds the symbol exponent s across the strip — '(C)=0 ⟹ WH index=0' is a CATEGORY ERROR"
    - "INHERITS the firewall: edge_bc=True drops the z_min node = enforces admissibility (B=0) by hand, so '0 unstable eigenvalues' restates the firewall, not independent evidence"
    - "CANNOT reach Re λ≥1: discrete Hilbert violates H[(Y-1)^{-1/2}]=0 at the edge node by −0.445·e^{−z_min/2} (reproduced to 3 sig figs, growing with depth) — out-of-strip uncertifiable"
  independent_rederivation: true
  machine_verification: "validate/edge_argprinciple_count.py (reproduced in workflow wf_e6371da0)"
lean_ref: none
tags: [track2-evidence, in-strip-only, NOT-the-index]
---

The cheap risk-bounding check (synthesis's first concrete step), done properly as a count.
Robust corroboration that no unstable point mode exists in the IN-STRIP region. Filed with
its three caveats front-and-center because conflating it with the Wiener–Hopf index would be
the NS-C6/SQG circular move (it is a different winding, inherits the firewall, and cannot reach
the decisive Re λ≥1 out-of-strip regime). Companion to the cleaner `ξ=(Y-1)^{1/2}` unit-recurrence
form of LEM-07 (`a_{m+1}=(m+λ)(m+2)tan(πm/2)a_m`, sympy-verified exact). The leap (LEM-06) stays
`reduced-to-sharp-WH-factorization`; this does NOT advance it. See `docs/desk/edge_mellin_symbol.md` §5.1.
