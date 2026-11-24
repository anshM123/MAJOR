#!/bin/bash
# EDGE_R2.md: reproduce (~1.3 CPU-h; at most 2 processes).  Read-only imports from edger/, highspin/, neargap2/, spinfunc/, agproof/,
# p4win/, lsw/; input ads4/ext/results/ext_results_batch3_partial.json (JSON only; no external pickle is loaded).
set -e
cd "$(dirname "$0")"
ulimit -v 3000000; export OMP_NUM_THREADS=1
printf '*\n!.gitignore\n' > cache/.gitignore
python3 e2_merge.py                                   # §1 batch-3 n1 rows -> cache/pts (8 states)
./run_k2.sh                                           # §1.2 K2 basis records + shadows (~45 min, one lane)
python3 e2_calib.py 26,30,34,38,42,46,50 > cache/calib.out &     # §3 exact-scalar calibration (~4 min)
for J in 26 30 34 33 38 42 46 50; do python3 e2_gamma.py $J > cache/gamma_J$J.out; done; wait   # §2 Gamma, costs (~1 min)
python3 e2_report.py                                  # table
