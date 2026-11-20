#!/bin/bash
# ee1000 queue A5 (after A4): the remaining odd-J K2 rows (35 nodes). Resumable.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 EE1000_Q=A
while pgrep -f ee1000_queueA4.sh > /dev/null; do sleep 30; done
for J in 29 25 21; do python3 ee1000_rows.py odd EE,C3H $J $((1000 + J)) X35; done
echo QUEUE_A5_DONE
