# Paper: Stability spectra of the singular Hou–Luo profile

**Status:** ready-to-arXiv (2026-07-22)  
**PDF:** `stability_spectra_singular_houluo.pdf`  
**Source:** `stability_spectra_singular_houluo.tex`  
**Legacy draft:** `stability_spectra_singular_houluo.md` (superseded by `.tex`/`.pdf`)

## Build PDF

```bash
cd paper
pdflatex stability_spectra_singular_houluo.tex
pdflatex stability_spectra_singular_houluo.tex   # TOC / refs
```

## Regenerate figures (one command each)

```bash
# from repo root
python3 paper/scripts/fig_spectrum.py
python3 paper/scripts/fig_gap.py
python3 paper/scripts/fig_coevolve.py
python3 paper/scripts/fig_nonlinear.py
```

Outputs: `paper/figures/fig_{spectrum,gap,coevolve,nonlinear}.{pdf,png}`

## Scoop log

`docs/scoop_recheck_2026_07_22_stability_spectra.md`
