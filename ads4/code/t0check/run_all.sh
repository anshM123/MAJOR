#!/bin/bash
# t0cert/run_all.sh (T0_CERT.md §6): full reproduction, <= 2 processes, ~ 20 CPU-min in total.
#   T0_DGMIN = smallest gap covered (the certificate holds for every Delta_g >= T0_DGMIN), T0_C = the constant c,
#   T0_BPOS = split point (default 16).  t3_flat's optional argument = the threshold actually used (>= t1's T_abs; the
#   slack absorbs the realised-weight eps_W term of T0_CERT §4.2).
set -e
cd "$(dirname "$0")"
ulimit -v 3000000; export OMP_NUM_THREADS=1
run () {   # DG C TFLAT TAG
  export T0_DGMIN=$1 T0_C=$2
  python3 t2_lrinf.py > cache/t2_$4.log 2>&1 &              # Part L  (bhat > 2 Lam)
  python3 t1_curv.py  > cache/t1_$4.log 2>&1                # Parts K, T (bhat <= 2 Lam; writes the threshold json)
  wait
  python3 t3_flat.py $3 > cache/t3_$4.log 2>&1              # Part F (bhat <= B_POS)
  grep -h VERDICT cache/t1_$4.log cache/t2_$4.log cache/t3_$4.log
}
run 1000 199/1000   3/100000    1000
run 300  199/1000   331/1000000 300
run 1000 1999/10000 ""          1000_c1999
python3 t0_checks.py > cache/t0_checks.log 2>&1             # [DV] sanity checks (not part of the proof)
