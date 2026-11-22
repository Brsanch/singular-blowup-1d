#!/bin/bash
cd "/Volumes/4TB SD/ClaudeCode/singular-blowup-1d"
export PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
L=/private/tmp/hdiag_review/logs_ad4; : > $L/_rc.txt
# the two gates edited in the diff are run from their ad4 text
git show ad4a0e5:validate/edge_eprime_probe.py > /private/tmp/hdiag_review/old_ad4/edge_eprime_probe.py
git show ad4a0e5:validate/stability_edge_z_spectrum.py > /private/tmp/hdiag_review/old_ad4/stability_edge_z_spectrum.py
for s in cap_scaffold energy_estimate_search edge_layer_obstruction spectral_gap theta_munu_gap edge_mellin_symbol edge_argprinciple_count edge_spectrum_continuation edge_w3_conjugation_check edge_coupled_resolvent_scan edge_gauge_a_eigenmode_check edge_nonlinear_bootstrap_closure edge_nonlinear_comoving_closure edge_coupled_h_nonlinear_closure adversarial_review_2026_11_21_probe coevolve_march nonlinear_shape_march coevolve_moving_frame coevolve_stress coevolve_paper_match coevolve_projected; do
  /Users/bryan/.claude/bin/reslock run --cpu 4 --purpose "hilbert diag review" -- python3 /private/tmp/hdiag_review/run_old.py "$PWD/validate/$s.py" > $L/$s.log 2>&1
  echo "$s rc=$?" >> $L/_rc.txt
done
echo DONE >> $L/_rc.txt
