---
id: LEM-06
title: (★) The nonlocal forcing is dominated by the anisotropic transport dissipation — ⟨G, ½Q_op G⟩_W ≤ ∫φ'WG² − δ‖G‖²_W for G ⊥_W 𝒩 with G(−∞)=0, for some anisotropic weight W=e^{-φ} (φ'→2γ in the bulk) and δ>0.
category: LEM
status: conjectured
confidence: 0.60
dependencies: [LEM-04, LEM-07, LEM-09, OBS-01, OBS-03, OBS-04, OBS-07]
supports: [THR-01]
cross_citations:
  - {label: "Chen 2026", ref: "arXiv:2605.15130", note: "finite-codimension + anisotropic + IBP-twice template"}
verification:
  method: numerics
  rounds: 3
  adversarial_vectors:
    - "single exponential weight φ=2γz FAILS (field-of-values diverges, C_force 10→83) — refuted; the anisotropic weight is required"
    - "BULK abscissa < 0 and depth-stable (−0.41/−0.65 across z_min=−10..−16) — the bulk genuinely closes"
    - "6-lens derive + 18 adversarial verifies (2026-06-13): the EDGE layer does NOT close by ANY weighted-energy/Schur/flux method — the forcing inherits the (Y-1)^{-1/2} profile amplitude at the edge (A_full·e^{z_min/2}=0.099 const, reproduced); the +40 lives on the g-unbounded branch (κ=1e-3 penalty crushes it to negative)"
    - "circularity firewall: H(g)(1)=0 is AUTOMATIC for L²_W modes (not a new lemma); asserting it from the benign spectrum = the NS-C6/SQG circular move — declined"
  independent_rederivation: true
  machine_verification: "validate/energy_estimate_search.py (bulk<0) + validate/edge_layer_obstruction.py (self-similar collapse 0.099·e^{−z_min/2}; depth-stable bulk; g-penalty crush; discrete-Hilbert −0.445·e^{−z_min/2}); NOT a closed analytic proof"
lean_ref: none
tags: [track2, leap, critical-submission-blocker]
---

**The non-interpolative leap — RESCOPED 2026-06-13 (`docs/desk/edge_layer_estimate.md`).**
(★) is numerically established in the BULK (depth-stable coercive δ=0.02..0.78). The
edge-layer leg was attacked by 6 independent lenses + 18 adversarial verifications and
the verdict is a sharp NEGATIVE result: **no bounded-at-edge weighted-energy / Schur /
flux method closes the edge layer** — the nonlocal forcing inherits the base-profile
`(Y-1)^{-1/2}` amplitude there, and the only weight that would dominate it (φ'≳1 at the
edge) re-triggers the convective field-of-values divergence (cap_scaffold). The `+40`
edge abscissa is a continuum non-normal profile-resonance carried by the g-unbounded
`(Y-1)^{-1/2}` indicial branch (`Re λ∈[−1,0)`: L²_W-admissible but `g=e^{−z/2}G∉L²(dz)`),
excluded by admissibility — but only enforced numerically, not proven.

So the leap is NO LONGER an energy estimate: it is the **non-self-adjoint spectral /
Fredholm** statement "the L²_W realization of `M_GG ⊖ 𝒩` has no `Re λ ≥ 0` spectrum".

**SHARPENED 2026-06-14 (LEM-07, `docs/desk/edge_mellin_symbol.md`).** The edge Mellin
symbol is now computed and 4-method verified: the half-line Hilbert is the multiplier
`𝓗(p)=cot(πp)`, and `½Q_op` is a `−1/2` Mellin ladder `c(ζ)=−tan(πζ)/(2ζ+1)`. The
eigenvalue problem reduces to the explicit **Mellin difference equation**
`(2s+λ)a(s)=c(s+1/2)a(s+1/2)`, and the ENTIRE remaining obligation is a single
**Wiener–Hopf index**: show the winding/partial-index of the factorization (edge line →
far-field line) is `0` for all `Re λ≥0` mod `𝒩`. In-strip `Re λ∈[0,1)` shows no automatic
pole (encouraging); the decisive `Re λ≥1` needs the analytic continuation. Probably a small
research problem (graded-lattice `cot`-symbol factorization with the `𝒩`-quotient), NOT an
off-the-shelf lemma. Circularity firewall: `B=a(0)` must be COMPUTED from the WH product,
never inferred from the benign spectrum. Stays `conjectured` (the index is named, not
evaluated).

**SHARPENED again 2026-06-14 (`docs/desk/edge_mellin_symbol.md` §5.1–5.2, OBS-05/06).** The
in-strip argument-principle count = 0 (OBS-05) is robust corroboration but NOT the index (a
different winding; inherits the firewall; can't reach `Re λ≥1`). The continuum Mellin-line run
(OBS-06) established: (i) lattice/power-series solutions are edge-divergent for ALL λ, so the
eigenfunctions are Mellin-line (`Ĝ_p(s)=Γ(s+λ)Γ(s+2)τ(s)`, `τ`=tan-Barnes); (ii) the symbol
winding `W=+3/2` (`Re λ>0`)/`+½` (`Re λ<0`) is the indicial-root TAUTOLOGY, not the Fredholm
index — the symbol VANISHES on the Plancherel line `Re s=0` (`R(0)=0 ∀λ` ⟹ non-Fredholm there),
and even index 0 ⊅ empty kernel (**Coburn's lemma blocked** by the on-line zero). The single
remaining object is now sharp: **construct the tan-Barnes `τ(s)` holomorphic + `L²`-decaying on
the correct (side-resolved) inversion contour, plus a Coburn-substitute / Pincherle-non-coincidence
step** ⟹ `dim ker = 0` for `Re λ≥0` mod `𝒩`. The genuine spectral content lives entirely in `τ(s)`.

**τ CONSTRUCTED — ELEMENTARY (2026-06-14, `edge_mellin_symbol.md` §5.3, machine-verified).** The
tan-Barnes is `τ(s)=π/sin(πs/2)=Γ(s/2)Γ(1−s/2)`, so `Ĝ_p(s)=Γ(s+λ)Γ(s+2)Γ(s/2)Γ(1−s/2)` (four
ordinary Gammas; eigenfunction = Fox‑H), with the WH split `Γ(s/2)`[left] `| Γ(1−s/2)`[right] read
straight off. This discharges the "construct τ" piece. The point spectrum reduces to an explicit
pole-killing-vs-`L²`-decay BALANCE, now made RIGOROUS (LEM-08, §5.4): the exact integer obstruction
`d=0` (line-`L²`) ∧ `d≥1` (kill the `s=0` edge-const) is empty ⟹ **no `Re λ≥0` POINT eigenvalue mod
`𝒩`** (machine-gated, firewall-clean, 3 routes agree). 

**So the POINT-spectrum half of (★) is resolved (LEM-08); LEM-06 is now the ESSENTIAL-spectrum half.**
The genuine remaining object: the essential spectrum sits on `Re s=0` (`R(0)=0 ∀λ` ⟹ non-Fredholm
Plancherel line), so a **strict** stability margin `Re σ≤−δ` (asymptotic stability, Conjecture 2.4)
still requires moving it off the imaginary axis — which needs the edge-*unbounded* convective weight
that re-triggers `cap_scaffold`'s field-of-values divergence. That is the original edge-weight wall,
now precisely isolated as the *only* thing between point-spectral stability and a strict gap. Secondary:
a clean Lean/interval statement of "bounded-at-edge `W` ⟹ const ∉ `L²_W`" (the computed firewall content).

**RESCOPED AGAIN 2026-07-15 (LEM-09, `docs/desk/convective_weight_spectrum.md`) — the
edge-weight wall is NOT a wall for the *spectral* problem.** Working directly in the convective
space `L²_σ = L²(e^{−2σz}dz)`, `σ∈(0,½)`: (i) the ε-similarity makes `M` unitarily scale-free in
the forcing, so the `cap_scaffold` FoV divergence was a truncation artifact *as a spectral
diagnostic* — while the energy route stays dead for a now-PROVEN reason (the exact `w=e^{−z/4}u`
diagonalization has a positive multiplier band at `σ+¼`, realized by the resonant tail
`e^{(σ+¼)z}`); (ii) the far-field transport line `{Re λ = −2σ} ⊆ σ_ess` caps the gap at `2σ`
(sup `1⁻` over the family); (iii) point spectrum right of the line is EMPTY (mod the named
edge-bootstrap Lemma E′), `𝒩` excluded by the space; (iv) the resolvent is EXPLICIT
(`Π_n = X/[X(·+n/2)(λ+2s+n)]`, defect-closed, index-0, machine-verified end-to-end). **The
essential-spectrum half of LEM-06 is therefore now ONE named estimate: the uniform `L²_σ` bound
of the defect-closed ladder resolvent on `Re λ ≥ −2σ+δ`** (≡ weak-coupling uniformity; probed by
OBS-07). Named route: bent-contour Barnes resummation + Schur test on the explicit Gamma-ratio
kernel (vertical contour provably diverges at rate `π/4` — the bend is mandatory). After that:
semigroup generation + Gearhart–Prüss on the constrained domain, then the `M_Gh` coupling, then
the nonlinear leg (all still open).
