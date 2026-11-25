#!/bin/bash
# t0cert/run_all.sh (T0_CERT.md): full reproduction, <= 2 processes, ~ 25 CPU-min in total.
#   T0_DGMIN = smallest gap covered (every Delta_g >= T0_DGMIN), T0_C = the constant c, T0_BPOS = split point (default 16).
set -e
cd "$(dirname "$0")"
ulimit -v 3000000; export OMP_NUM_THREADS=1
for DG in 1000 300; do
  export T0_DGMIN=$DG T0_C=199/1000
  python3 t2_lrinf.py  > cache/t2_$DG.log 2>&1 &          # Part L  (bhat > 2 Lam)
  python3 t1_curv.py   > cache/t1_$DG.log 2>&1            # Parts K, T (bhat <= 2 Lam; writes the threshold json)
  wait
  python3 t3_flat.py   > cache/t3_$DG.log 2>&1            # Part F (bhat <= B_POS; reads the threshold from t1's json)
  grep VERDICT cache/t1_$DG.log cache/t2_$DG.log cache/t3_$DG.log
done
python3 t0_checks.py > cache/t0_checks.log 2>&1           # [DV] sanity checks (not part of the proof)
