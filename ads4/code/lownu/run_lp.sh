#!/bin/bash
# LF1/LF2 scan: LP + dense check for a in {0.025,0.05,0.075,0.1}, delta in {0.2, 0.1}, matter rows off/on.  ~15 min, 1 process.
ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$(dirname "$0")"
for a in 0.025 0.05 0.075 0.1; do for m in 0 1; do for dl in 0.2 0.1; do
  python3 ln2_lp.py $a 5 8 $dl $m | tail -1 | tee -a cache/run_lp.txt
done; done; done
for f in cache/lp_a*_dl0.20_m*.json; do python3 ln3_check.py $f | tee -a cache/run_check.txt; done
