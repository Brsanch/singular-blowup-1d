---
id: LEM-09
title: Convective-weight structure pack for M_GG = −2∂_z + ½Q_op on L²_σ = L²(e^{−2σz}dz), σ∈(0,½) — (A) ε-similarity (M unitarily ≅ −2∂_z + εK ∀ε>0; spectrum/resolvent norms/FoV scale-free, so truncated FoV divergences are similarity-breaking artifacts); (B) far-field transport line {Re λ = −2σ} ⊆ σ_ess, capping the gap at 2σ (sup 1⁻); (C) energy no-go with exact mechanism: w = e^{−z/4}u diagonalizes the K-form to the multiplier Re c(σ+¼+it), which has a positive band (Re c(σ+¼)>0 ∀σ∈(0,½)) realized by the resonant tail e^{(σ+¼)z} — no pure-exponential weight is dissipative; (D) rank-one defect: KG = (L²_σ-part) − (𝒢(½)/2π)𝟙, so D(M) ⊆ {𝒢(½)=0} (the B=H(g)(1)=0 firewall IS the generator's domain condition), M(X-inv) = λ(X-inv) − Γ(1+λ)𝟙, and Res₀𝒢 = 𝒢(½)/(2πλ) (domain condition ⟺ no far-field constant); (E, conditional on E′) no L²_σ point spectrum in Re λ > −2σ, 𝒩 included (excluded by the space); (F) exact resolvent: Π_n(s) = X(s)/[X(s+n/2)(λ+2s+n)], defect-closed ladder solves (λ−M)G = F on the dense entire-transform class, one periodic freedom ↔ one defect condition (index-0).
category: LEM
status: proved-on-paper
confidence: 0.8
dependencies: [LEM-07, LEM-08]
supports: [LEM-06, THR-01]
cross_citations:
  - {label: "X in z-form", ref: "X(s)=Γ(2s+λ)Γ(2s+2)Γ(s)Γ(1−s), (λ+2s)X(s)=c(s+½)X(s+½)", note: "LEM-08's 4-Gamma at doubled argument; note s_ξ↔s_z conventions REFLECT — end assignments re-derived, not imported"}
  - {label: "variation of parameters", ref: "Π_n(s)=X(s)/[X(s+n/2)(λ+2s+n)]", note: "closed form; makes the resolvent a Barnes-type object"}
verification:
  method: proof
  rounds: 1
  adversarial_vectors:
    - "functional equation machine-exact (2e-31, complex λ,s); ladder symbol K e^{sz}=c(s)e^{(s−½)z} re-verified by pure continuum quadrature (9e-14) — no cot-table shortcut"
    - "defect SIGN measured, not assumed: K(bump) → −𝒢(½)/(2π) at the edge (−0.182141 vs −0.182134); the naive-eigenfunction scenario (spectrum fills the half-plane) explicitly killed by the machine-verified defect identity MG_p = λG_p − Γ(1+λ)𝟙 (rel 3e-4; far-const Γ(λ); edge decay exactly e^z)"
    - "resolvent verified END-TO-END: line equation 5e-15, operator PDE residual 1.4e-3 (quadrature-limited), Ĝ(½)=0 by construction, far-field constant 4e-5 ≈ 0 (residue⇔value identity X(½)/(2πλ)=Γ(λ)=Res₀X exact)"
    - "s=½ is a removable zero×pole point of the ladder (odd n vanish, even n have finite c(1)c(3/2)-type limits, −1/12 at n=2): naive 𝒢₀(½) evaluation is WRONG (14.25 vs true 4.30 in the test case) — cost one bug, now limit-evaluated"
    - "energy no-go artifact-free: Re c(σ+¼)>0 (7.79/0.85/0.43 at 0.52/0.65/0.74), resonant-tail form +0.18/+1.09/+4.04 at L=4/8/16 (divergent positive); localized bumps DISSIPATIVE (−0.579, exact e^{−a/2} scaling to 8 digits) — the cap_scaffold C_force divergence is superseded (it conflated this continuum σ+¼ resonance with the discrete-Hilbert edge bias, same e^{−z_min/2} rate)"
  independent_rederivation: false
  machine_verification: "validate/edge_convective_structure.py (12 checks, GATE PASS 2026-07-15)"
lean_ref: none
tags: [track2, convective-weight, resolvent, LEM-06-attack]
---

**CORR 2026-11-21 (LEM-10):** part (C) holds only for `σ > ¼`. `Re c(σ+¼) > 0 ∀σ∈(0,½)` in the title is
wrong (`c < 0` on `(0,½)`; the sampled values all have `σ ≥ 0.27`). For `σ ≤ ¼` the `L²_σ` energy estimate
closes with rate `2σ` (manuscript Proposition `prop:energy_form`). The title is left as written for the record.
Also superseded: "spectrum = the line" (it is the closed half-plane `Re λ ≤ −2σ`, LEM-10).

Full derivations + proofs: `docs/desk/convective_weight_spectrum.md` (2026-07-15). One honest
conditional: **(E) rests on Lemma E′** (edge bootstrap: an `L²_σ` eigenfunction has no
`e^{−λz/2}` homogeneous branch, hence is strip-regular — the continuum form of the repo's
inflow-BC/admissibility step, named and open). Everything else is identity-grade
(machine-verified 13–31 digits) or short-proof-grade (ε-similarity, Liouville transfer with the
¾-coefficient bound, Weyl sequences).

**What it buys for LEM-06:** the essential-spectrum half is reduced to ONE estimate — a uniform
`L²_σ` bound on the defect-closed ladder resolvent for `Re λ ≥ −2σ+δ` (≡ weak-coupling
uniformity by (A); probed artifact-free by OBS-07). The named analytic route is the bent-contour
Barnes resummation of `Σ_n 𝓕(s+n/2)/[X(s+n/2)(λ+2s+n)]` + Schur test (vertical-line rotation
provably diverges at rate π/4 — measured 0.73 → π/4; the bend is mandatory). NOT claimed: the
bound itself, semigroup generation, Gearhart–Prüss, the `M_Gh` coupling, anything nonlinear.
