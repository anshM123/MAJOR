#!/bin/bash
# ee1000 queue A2 (replaces the rest of A): even-J K2 rows at Delta = 1000 + J with 41 x nodes (0.025 ... 0.9). Resumable.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 EE1000_Q=A
while pgrep -f "ee1000_rows.py even EE,C3H 40 1040 X39" > /dev/null; do sleep 20; done
python3 ee1000_rows.py even EE,C3H 40 1040 X41
python3 ee1000_rows.py refb 48,56,64
python3 ee1000_rows.py even EE,C3H 32 1040 LOHI
python3 ee1000_rows.py even EE,C3H 36 1036 X41
python3 ee1000_rows.py even EE,C3H 48 1048 X41
python3 ee1000_rows.py even EE,C3H 64 1064 X41
python3 ee1000_rows.py even EE,C3H 56 1056 X41
echo QUEUE_A2_DONE
