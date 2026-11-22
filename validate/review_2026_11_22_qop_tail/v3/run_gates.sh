#!/bin/bash
cd "/Volumes/4TB SD/ClaudeCode/singular-blowup-1d"
export PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
: > /private/tmp/hdiag_review/logs/_rc.txt
for s in stability_edge_z_spectrum edge_eprime_probe cap_scaffold energy_estimate_search edge_layer_obstruction spectral_gap theta_munu_gap edge_mellin_symbol edge_argprinciple_count edge_spectrum_continuation edge_w3_conjugation_check edge_coupled_resolvent_scan edge_gauge_a_eigenmode_check edge_nonlinear_bootstrap_closure edge_nonlinear_comoving_closure edge_coupled_h_nonlinear_closure adversarial_review_2026_11_21_probe coevolve_march nonlinear_shape_march coevolve_moving_frame coevolve_stress coevolve_paper_match coevolve_projected; do
  t0=$(date +%s)
  /Users/bryan/.claude/bin/reslock run --cpu 4 --purpose "hilbert diag review" -- python3 validate/$s.py > /private/tmp/hdiag_review/logs/$s.log 2>&1
  rc=$?
  echo "$s rc=$rc secs=$(( $(date +%s)-t0 ))" >> /private/tmp/hdiag_review/logs/_rc.txt
done
echo DONE >> /private/tmp/hdiag_review/logs/_rc.txt
