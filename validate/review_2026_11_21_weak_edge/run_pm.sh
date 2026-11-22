#!/bin/bash
export PYTHONDONTWRITEBYTECODE=1
cd "$(dirname "$0")"
run() { /Users/bryan/.claude/bin/reslock run --cpu 4 --purpose "weak-edge review probe" -- python3 pm_analysis.py "$@" > "pm_$1_$2_$3_$4_$5_$6.log" 2>&1; echo "=== $*"; grep -v "^#   c_adv" "pm_$1_$2_$3_$4_$5_$6.log"; }
run 0.6 -48 48 4916 u3 0.02
run 0.6 -104 24 6554 u3 0.02
run 0 -104 24 6554 u3 0.02
run 0.6 -104 24 13108 u3 0.02
