#!/bin/bash
# BELOW_FOLD2.md: production chain (one core each, <= 3 GB).  bg_eeq ~2 min (J=0) / ~14 min (J=2) per state incl. reference.
cd "$(dirname "$0")"
ulimit -v 3000000
export OMP_NUM_THREADS=1
python3 bg_scalar.py 6 && python3 bg_validate.py 6                        # scalar K = 6 moments + validation table
python3 bg_scalarJ.py 14 2 && python3 bg_scalarJ.py 12 4 && python3 bg_validateJ.py && python3 bg_lrscalar.py
for st in "300 0" "1000 0" "314.925 0" "314.935 2" "1049.93 2"; do
  set -- $st
  [ -f cache/bg_eeq_D$1_J$2_K4.pkl ] || python3 bg_eeq.py $1 $2 4 5 2 > cache/bg_eeq_D$1_J$2.log 2>&1
done
python3 bg_lrcheck.py 300 0 4; python3 bg_eetrunc.py 300 0 4; python3 bg_eetrunc.py 1000 0 4; python3 bg_eetrunc.py 314.935 2 4
python3 bg_qan.py 300 0 4 300 ++; python3 bg_qan.py 1000 0 4 1000 ++; python3 bg_qan.py 314.925 0 4 300 ++
python3 bg_qan.py 314.935 2 4 300 ++; python3 bg_qan.py 1049.93 2 4 1000 ++
python3 bg_covfit.py 300 0; python3 bg_covfit.py 1000 0; python3 bg_covfit.py 314.935 2; python3 bg_ffscan.py
python3 bg_k2check.py 0      # K2 Laplace rows (w1_pts) at nu = 19.9, 39.9: inconclusive (Laplace unconverged), BELOW_FOLD2 §3.3
