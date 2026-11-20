#!/bin/bash
# ee1000 queue A6 (after A4): a K2 row for J = 64 at a mass where K2 converges (Delta = 2900, m = 2.9 Lam, J^2/m = 1.41), nodes x = 0.05 ...
# 0.35 (step 0.0125; the states there have x <= 1/2.9): measures the beyond-Regge correction where the band model is used for J = 64.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 EE1000_Q=A
while pgrep -f ee1000_queueA4.sh > /dev/null; do sleep 30; done
python3 ee1000_rows.py even EE,C3H 64 2900 $(python3 -c "print(','.join('%g' % (0.05 + 0.0125*i) for i in range(25)))")
echo QUEUE_A6_DONE
