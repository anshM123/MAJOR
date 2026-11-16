#!/bin/sh
# Regenerates everything for ads4/SKELETON_TEST.md.  Uses at most 2 cores (two chains in parallel).
# Imports, unchanged, ads4/code/e2e, ads4/code/fourpt (incl. fourpt/tttt/t2_mft_heavy_l1.factor), ads4/code/uniform (CHMRS series:
# needs uniform/cache/chmrs_series.pkl or the ancillary file, see UNIFORM_REMAINDER.md).  Caches/outputs: cache/ (git-ignored).
# Block tensors: ~1.5 min (m ~ 10-17) ... ~35 min (m ~ 70) per case on one core; the fourpt cache (cache/modes4 there) is reused when present.
# 'sh run_all.sh fast' skips the block-tensor computations (uses whatever is cached) and only redoes the analysis.
set -e
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1   # <= 2 cores in total
cd "$(dirname "$0")"
mkdir -p cache/logs
if [ "$1" != "fast" ]; then
  # precision: scalar moments on refined grids; two-pair moments on the u-refined grid at m ~ 10-12 (agreement with the standard grid)
  python3 s0_scalar_precision.py 14.166224 2 2.0 > cache/logs/s0.log
  python3 s1_spin_scan.py 14.166224 2 u2 12 1 > cache/logs/s1_D14_u2.log &
  python3 s1_spin_scan.py 13.401744 5 u2 12 1 > cache/logs/s1_D13_u2.log &
  wait
  # two chains (beta ~ 0.4 and beta ~ 1), standard grid, nu = x m, x = 0 ... 0.65 in steps of 0.05 (nu <= 16)
  ( for c in "14.166224 2" "19.232713 3" "24.299203 4" "34.432182 6" "44.565161 8" "74.964099 14" "105.363036 20"; do
      python3 s1_spin_scan.py $c > "cache/logs/s1_$(echo $c | tr ' ' '_').log"; done ) &
  ( for c in "13.401744 5" "19.893604 8" "24.221511 10" "35.041278 15" "45.861045 20" "80.484300 36" "115.107554 52"; do
      python3 s1_spin_scan.py $c > "cache/logs/s1_$(echo $c | tr ' ' '_').log"; done ) &
  wait
fi
python3 s4_summary.py > cache/out_s4_summary.txt
