#!/bin/bash
# Usage: ./reproduce.sh verify   (≈3 min: re-runs 1 seed of every Phase-3 block and re-runs all analysis scripts, diffs against saved results)
#        ./reproduce.sh full     (hours: regenerates every Phase-3 result pickle from scratch; see README_REPRODUCE.md for times)
set -e; cd "$(dirname "$0")"
case "$1" in
verify)
  python3 repro_check.py s25 s26 s27 s28 s30            # expect: max abs diff 0.000e+00 for every block
  for a in s25 s26 s27 s28 s30; do python3 phase3_${a}_analyze.py > /tmp/re_$a.txt; diff -q /tmp/re_$a.txt phase3_${a}_analysis.txt && echo "$a analysis identical"; done
  python3 phase3_pilot_analyze.py stage1 | diff -q - phase3_pilot_stage1_analysis.txt && echo "stage1 analysis identical"
  python3 phase3_cdf_diag_analyze.py | diff -q - phase3_cdf_diag_analysis.txt && echo "cdf diag analysis identical";;
full)
  python3 phase3_pilot.py 1 30 stage1 && python3 phase3_pilot_analyze.py stage1      # ~70 min, 2 procs (original run resumed twice; resumable)
  python3 phase3_cdf_diag.py                                                          # ~9 min
  python3 phase3_s25_pilot.py && python3 phase3_s25_analyze.py                        # ~5 min
  python3 phase3_s26.py 31 60 phase3_s26.pkl && python3 phase3_s26_analyze.py         # ~2 min
  python3 phase3_s27.py 61 90 phase3_s27.pkl && python3 phase3_s27_analyze.py         # ~20 s
  python3 phase3_s28.py 301 330 phase3_s28.pkl && python3 phase3_s28_analyze.py       # ~1 s
  python3 phase3_s28.py 1001 1300 phase3_s30.pkl && python3 phase3_s30_analyze.py     # ~1 min  (TEST SEEDS: H1 only)
  python3 export_csv.py;;
*) echo "usage: $0 verify|full"; exit 1;;
esac
