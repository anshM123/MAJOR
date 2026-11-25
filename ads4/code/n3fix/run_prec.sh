#!/bin/bash
# N3C3H_FIX.md §3.3: precision check of the Gamma-level n3 / n1 difference at J = 26, x = 0.92 (n3 M = 12 with an independent seed,
# n1 M = 13).  Started after run_more.sh; two processes.
cd "$(dirname "$0")"
ulimit -v 3000000; export OMP_NUM_THREADS=1
while pgrep -f "nf_run.py (edge3|ovl) " > /dev/null; do sleep 20; done
python3 nf_run.py prec_n3 n3m12 26 1026 0.92 7000 > cache/prec_n3_J26.out 2>&1 &
python3 nf_run.py prec_n1 n1m13 26 1026 0.92 > cache/prec_n1_J26.out 2>&1
wait
