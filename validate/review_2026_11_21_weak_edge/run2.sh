#!/bin/bash
export PYTHONDONTWRITEBYTECODE=1
cd "$(dirname "$0")"
/Users/bryan/.claude/bin/reslock run --cpu 4 --purpose "weak-edge review probe" -- python3 pm_scan.py "$1" > "pm_scan_$1.log" 2>&1
cat "pm_scan_$1.log"
