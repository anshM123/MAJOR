#!/bin/bash
# Delta = 50: EE D' grids and n3 crystal runs (two seeds) at nu = 12, 14, 16; runs after farw_q.sh
ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$(dirname "$0")"
while pgrep -f "farw_q.sh" > /dev/null; do sleep 10; done
for nu in 16 12 14; do
  [ -f cache/cr_D50_nu${nu}.pkl ] || python3 pcyc_crys.py 50 $nu 12 8 > cache/log_cr_D50_nu${nu}.txt 2>&1
  python3 pcyc_cv.py cache/cr_D50_nu${nu}.pkl 2 >> cache/log_cr_D50_nu${nu}.txt 2>&1
done
for g in 12,12,6,12,12,20,12 12,14,8,12,14,24,16; do
  [ -f cache/dq_D50_E_${g//,/-}.npz ] || python3 pcyc_dquad.py 50 E $g > cache/log_dq_50_E_$g.txt 2>&1
done
for nu in 16 12 14; do
  [ -f cache/cr_D50_nu${nu}_s2.pkl ] || PC_SEED=5151 PC_TAG=_s2 python3 pcyc_crys.py 50 $nu 12 8 > cache/log_cr_D50_nu${nu}_s2.txt 2>&1
  python3 pcyc_cv.py cache/cr_D50_nu${nu}_s2.pkl 2 >> cache/log_cr_D50_nu${nu}_s2.txt 2>&1
done
