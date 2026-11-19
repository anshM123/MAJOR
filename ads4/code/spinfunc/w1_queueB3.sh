#!/bin/bash
# w1 queue B3 (after B2): direct EE at (336, 18) for larger x (K2 validation at J^2/m ~ 1 beyond x = 0.2)
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
while pgrep -f w1_queueB2.sh > /dev/null; do sleep 20; done
python3 w3_direct.py run 18 336 0.25,0.3,0.4,0.5
