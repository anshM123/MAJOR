#!/bin/bash
# fw queue B4 (supersedes fw_queueB3.sh after its J = 48 rows): J = 48 refb in the patched basis, extra reference-basis points
# x_ref = 0.15 ... 0.45 (the order-0 reference action is asymmetric above x_ref ~ 0.5 at J = 100), then the remaining bases.  Resumable.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
while pgrep -f "^python3 fw_rows.py" > /dev/null; do sleep 20; done
python3 fw_rows.py row EE,C3H 48 1048 X35
python3 fw_rmrefb48.py
python3 fw_rows.py refb EE 48
for J in 100 48 64; do python3 fw_rows.py state EE $J 20000 0.15,0.2,0.25,0.35,0.45; done
python3 fw_rows.py state EE 64 1064 0.15,0.2,0.25,0.3,0.35,0.4
python3 fw_rows.py refb EE 64
python3 fw_rows.py state EE 100 1100 0.15,0.2,0.25,0.3,0.35,0.4
python3 fw_rows.py state EE 64 364 0.15,0.2,0.25,0.3,0.35,0.4
for J in 150 200; do python3 fw_rows.py state EE $J 20000 0.15,0.2,0.25,0.3,0.35,0.4,0.45; done
python3 fw_rows.py state EE 150 450 0.15,0.2,0.25,0.3,0.35,0.4
python3 fw_rows.py state EE 200 1200 0.15,0.2,0.25,0.3,0.35,0.4
echo FW_QUEUE_B4_DONE
