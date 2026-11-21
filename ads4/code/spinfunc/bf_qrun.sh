#!/bin/bash
# BELOW_FOLD.md production: EE k=4 / scalar state-ref sums on the real section D (2 processes max)
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
(python3 bf_q.py run 300 0 17 16; python3 bf_q.py run 314.93453064396084 2 17 16) > cache/bf_qrun_a.log 2>&1 &
(python3 bf_q.py run 1000 0 17 16; python3 bf_q.py run 1049.927849210426 2 17 16) > cache/bf_qrun_b.log 2>&1 &
wait
