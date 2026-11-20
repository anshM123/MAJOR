#!/bin/bash
# ee1000 queue B (EE_DG1000_LARGEJ.md): direct-RQMC validation points, odd-J K2 rows at Delta = 1000 + J (39 x nodes).
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 EE1000_Q=B
python3 ee1000_direct.py 40 1040 0.6,0.5,0.7
python3 ee1000_rows.py odd EE,C3H 17 1017 X39
python3 ee1000_rows.py odd EE,C3H 31 1031 X39
python3 ee1000_direct.py 64 1064 0.3,0.6,0.15
python3 ee1000_rows.py odd EE,C3H 23 1023 X39
python3 ee1000_rows.py odd EE,C3H 19 1019 X39
python3 ee1000_rows.py odd EE,C3H 27 1027 X39
python3 ee1000_rows.py odd EE,C3H 21 1021 X39
python3 ee1000_rows.py odd EE,C3H 25 1025 X39
python3 ee1000_rows.py odd EE,C3H 29 1029 X39
echo QUEUE_B_DONE
