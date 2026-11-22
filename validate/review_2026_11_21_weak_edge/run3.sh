#!/bin/bash
export PYTHONDONTWRITEBYTECODE=1
cd "$(dirname "$0")"
for w in "$@"; do
/Users/bryan/.claude/bin/reslock run --cpu 4 --purpose "weak-edge review probe" -- python3 fs3_sweep.py "$w" > "fs3_$w.log" 2>&1
echo "##### $w"; cat "fs3_$w.log"; done
