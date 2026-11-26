# Cryptographic Priority Timestamps — Target A (1D Hou–Luo) & Target B (gCLM)

This directory contains an offline-verifiable cryptographic priority record for both manuscripts, Lean 4 formalizations, and validation gates in `singular-blowup-1d`:

- **`CHECKSUMS_TARGET_B.sha256`** — SHA-256 digests of:
  - `paper/stability_spectra_singular_gclm.pdf` & `.tex` (Target B: gCLM, 21 pages, unconditional global existence via Appendix A)
  - `paper/stability_spectra_singular_houluo.pdf` & `.tex` (Target A: 1D Hou–Luo, 22 pages, untouched)
  - `lean/GclmSingularBackbone.lean` (Lean 4.29 / Mathlib, `0 sorry`, standard axioms only)
  - `validate/gclm_{verify_singular_steady,linear_stability,branch2_and_dns_probe,rigorous_closure_gate}.py` (Gates 1–4, all PASS)
  - `docs/proof_gclm_singular_stability.md` & `docs/findings_gclm_singular_stability.md`
- **`CHECKSUMS_TARGET_B.sha256.tsq`** — RFC 3161 SHA-256 timestamp query generated from `CHECKSUMS_TARGET_B.sha256`.
- **`CHECKSUMS_TARGET_B.freetsa.tsr`** (with `freetsa_tsa.crt` and `freetsa_cacert.pem`) — RFC 3161 Trusted Timestamping Authority response signed by FreeTSA.org on **2026-11-26 00:40:07 GMT**.
- **`CHECKSUMS_TARGET_B.digicert.tsr`** — Independent backup RFC 3161 timestamp response signed by DigiCert TSA.

## Offline Verification Commands

Run from the repository root:

```bash
# 1. Verify file integrity against the SHA-256 manifest
shasum -a 256 -c docs/timestamps/CHECKSUMS_TARGET_B.sha256

# 2. Verify the FreeTSA RFC 3161 cryptographic timestamp signature
openssl ts -verify \
  -data docs/timestamps/CHECKSUMS_TARGET_B.sha256 \
  -in docs/timestamps/CHECKSUMS_TARGET_B.freetsa.tsr \
  -CAfile docs/timestamps/freetsa_cacert.pem \
  -untrusted docs/timestamps/freetsa_tsa.crt

# 3. Inspect exact UTC timestamp embedded in either token
openssl ts -reply -in docs/timestamps/CHECKSUMS_TARGET_B.freetsa.tsr -text | grep -E "Time stamp:|Status:"
openssl ts -reply -in docs/timestamps/CHECKSUMS_TARGET_B.digicert.tsr -text | grep -E "Time stamp:|Status:"
```
