#!/bin/bash
# fw queue B3 (supersedes fw_queueB.sh): fixed-engine K2 records in priority order: J = 48 rows (1048, X35, EE + C3H, state + reference),
# bases for (1064, 64), (1100, 100), (364, 64), (450, 150), (1200, 200) and their references.  Resumable.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
while pgrep -f "^python3 fw_rows.py" > /dev/null; do sleep 20; done
python3 fw_rows.py refb EE 100
python3 fw_rows.py row EE,C3H 48 1048 X35
python3 fw_rows.py state EE 64 1064 0.15,0.2,0.25,0.3,0.35,0.4
python3 fw_rows.py refb EE 64
python3 fw_rows.py state EE 100 1100 0.15,0.2,0.25,0.3,0.35,0.4
python3 fw_rows.py state EE 64 364 0.15,0.2,0.25,0.3,0.35,0.4
python3 fw_rows.py state EE 150 450 0.15,0.2,0.25,0.3,0.35,0.4
python3 fw_rows.py refb EE 150
python3 fw_rows.py state EE 200 1200 0.15,0.2,0.25,0.3,0.35,0.4
python3 fw_rows.py refb EE 200
echo FW_QUEUE_B3_DONE
