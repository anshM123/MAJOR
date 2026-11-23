#!/bin/bash
# crystal-cycle runs (pcyc_crys) + control-variate estimates (pcyc_cv); one chain.  Waits for the Delta = 30 D' chain to finish.
ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$(dirname "$0")"
while pgrep -f "dq_prod.sh 30" > /dev/null; do sleep 20; done
# above the fold: n3 hybrid contour (nu = 10 forced to n3 with PC_NUSPLIT = 9.9 as a cross-check of the n4 thimble sum)
for D in 30 100; do
  for nu in 12 14 16; do
    [ -f cache/cr_D${D}_nu${nu}.pkl ] || python3 pcyc_crys.py $D $nu 11 8 > cache/log_cr_D${D}_nu${nu}.txt 2>&1
    python3 pcyc_cv.py cache/cr_D${D}_nu${nu}.pkl 2 >> cache/log_cr_D${D}_nu${nu}.txt 2>&1
  done
  [ -f cache/cr_D${D}_nu10_n3.pkl ] || PC_NUSPLIT=9.9 PC_TAG=_n3 python3 pcyc_crys.py $D 10 11 8 > cache/log_cr_D${D}_nu10_n3.txt 2>&1
  python3 pcyc_cv.py cache/cr_D${D}_nu10_n3.pkl 2 >> cache/log_cr_D${D}_nu10_n3.txt 2>&1
done
# below the fold: n4 A+ + A- thimbles
for job in "30 10" "30 8" "30 4" "100 8" "30 3"; do
  set -- $job
  [ -f cache/cr_D$1_nu$2.pkl ] || python3 pcyc_crys.py $1 $2 10 8 > cache/log_cr_D$1_nu$2.txt 2>&1
  python3 pcyc_cv.py cache/cr_D$1_nu$2.pkl 2 >> cache/log_cr_D$1_nu$2.txt 2>&1
done
