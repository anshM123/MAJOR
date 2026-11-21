#!/bin/sh
# runs after stage3_run.sh: smeared-part exact check (Dg=100) and the dense off-grid check of the delta = 0.1/0.2 LP functionals
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
while pgrep -f stage3_run.sh > /dev/null; do sleep 20; done
python3 stage3_smeared_check.py > cache/stage3_smeared_check.log 2>&1 &
n=$(grep -c . cache/lp_g4.jsonl)
python3 dense_check.py $((n-3)) $((n-2)) > cache/stage3_dense.log 2>&1
wait
