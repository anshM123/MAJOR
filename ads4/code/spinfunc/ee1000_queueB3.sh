#!/bin/bash
# ee1000 queue B3 (replaces the rest of B2): direct-RQMC validation points (state; K2 reference), odd-J K2 rows at Delta = 1000 + J (35 nodes).
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 EE1000_Q=B
while pgrep -f "ee1000_direct.py 40 1040" > /dev/null; do sleep 20; done
python3 ee1000_rows.py odd EE,C3H 17 1017 X35
python3 ee1000_rows.py odd EE,C3H 31 1031 X35
python3 ee1000_direct.py 64 1064 0.3,0.6,0.15
python3 ee1000_rows.py odd EE,C3H 23 1023 X35
python3 ee1000_rows.py odd EE,C3H 27 1027 X35
python3 ee1000_rows.py odd EE,C3H 19 1019 X35
python3 ee1000_rows.py odd EE,C3H 21 1021 X35
python3 ee1000_rows.py odd EE,C3H 25 1025 X35
python3 ee1000_rows.py odd EE,C3H 29 1029 X35
echo QUEUE_B3_DONE
