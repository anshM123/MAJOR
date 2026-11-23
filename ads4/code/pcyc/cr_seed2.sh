#!/bin/bash
# independent-seed n3 runs at Delta = 100 (2x points, seed 5151); runs after exact_more.sh
ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$(dirname "$0")"
while pgrep -f "exact_more.sh" > /dev/null; do sleep 10; done
for nu in 16 14 12; do
  [ -f cache/cr_D100_nu${nu}_s2.pkl ] || PC_SEED=5151 PC_TAG=_s2 python3 pcyc_crys.py 100 $nu 12 8 > cache/log_cr_D100_nu${nu}_s2.txt 2>&1
  python3 pcyc_cv.py cache/cr_D100_nu${nu}_s2.pkl 2 >> cache/log_cr_D100_nu${nu}_s2.txt 2>&1
done
[ -f cache/cr_D100_nu10_n3s2.pkl ] || PC_SEED=5151 PC_NUSPLIT=9.9 PC_TAG=_n3s2 python3 pcyc_crys.py 100 10 12 8 > cache/log_cr_D100_nu10_n3s2.txt 2>&1
python3 pcyc_cv.py cache/cr_D100_nu10_n3s2.pkl 2 >> cache/log_cr_D100_nu10_n3s2.txt 2>&1
