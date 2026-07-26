# W2 — Lemma E′ (edge bootstrap / uniqueness)

**Status:** PARKED AS NAMED HYPOTHESIS — not falsified; not proved  
**Date:** 2026-07-25 (magnet #4 Phase 2)  
**Does not claim:** E′ proved · Prop E unconditional · INTERNAL-CAP-CERTIFIED

---

## Statement (desk §5 / Prop E)

**Lemma E′.** An \(L^2_\sigma\) eigenfunction of \(M\) for
\(\operatorname{Re}\lambda>-2\sigma\) has **no** homogeneous edge branch
\(G\sim e^{-\lambda z/2}\) (that weight diverges at the edge); it carries
recessive \(e^{z}\) edge decay and is strip-regular on
\([\sigma,\sigma+\tfrac12]\).

Prop E (no point spectrum right of the transport line) is **modulo E′**.

---

## Pre-registered claim B2

> E′ admits an analytic or CAP-backed argument consistent with the residual-floor
> gate (no constructible approx eigenfunction).

**Outcome:** **not cashed.** Residual floor is consistent with E′
(NUMERIC-FORECAST) but is not a continuum bootstrap proof. No circular
assumption introduced. Wall remains the named hypothesis E′.

---

## NUMERIC-FORECAST evidence (not a proof)

Gate `validate/edge_eprime_probe.py` · findings
`docs/findings_session_2026_07_22_eprime_residual_floor.md`.

At \(\sigma=0.4\), \(\lambda=0.3\):
```
KEY gamma_w=1.10035 r_min=0.585763 lam_spur_re=1.7066
    alpha_free=-0.544389 alpha_bc=3.02852 maxRe_bc=-0.00156221
```

| Fact | Reading |
|---|---|
| \(r_{\min}=0.586\); \(r(-14)/r(-6)=55.6\) | No constructible approx eigenfunction from \(G_{\mathrm{hom}}\); residual **grows** with depth |
| \(\gamma_w=\operatorname{Re}\lambda+2\sigma\) | Homogeneous branch not in \(L^2_\sigma\) (measured) |
| Free-edge spur \(\operatorname{Re}\lambda\approx +1.71\) killed by inflow BC | Discrete BC matches the continuum object E′ names |

**Falsification probe:** a continuum counterexample to E′ would need a sequence
with residual \(\to 0\) as edge depth \(\to\infty\). Observed discrete residual
moves the wrong way — supports keeping E′, does not prove it.

---

## Why not CAP this session

E′ is a boundary-regularity / Liouville-type statement. No fixed-point
formulation with computable defect bounds is ready. Opening CAP without
conditioning would violate “condition before certifying.”

---

## Park decision

- Keep E′ as a **named hypothesis** in Prop E.
- Do not weaken Prop E to hide the gap.
- Do not cite residual-floor KEY as PROVED-ON-PAPER.
- Reopen W2 only with a continuum Mellin/Liouville outline or a CAP-ready
  formulation — not more discrete residual scans.
