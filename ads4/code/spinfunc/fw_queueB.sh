#!/bin/bash
# fw queue B (FOLD_WINDOW_LARGEJ.md): fixed-engine K2 basis records (EE, state x = 0.15..0.8 and reference refb) for the large-J rows,
# then the J = 48 rows at (1048, 48) recomputed (EE + C3H, X35, state + reference).  Needs cache/fw_shadow_k2.pkl first (fw_shadow.py).
# Resumable (finished records skipped).
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 fw_rows.py state EE 100 400 B
python3 fw_rows.py refb EE 100
python3 fw_rows.py state EE 150 450 B
python3 fw_rows.py refb EE 150
python3 fw_rows.py row EE,C3H 48 1048 X35
python3 fw_rows.py state EE 64 364 B
python3 fw_rows.py state EE 64 1064 B
python3 fw_rows.py refb EE 64
python3 fw_rows.py state EE 100 1100 B
python3 fw_rows.py state EE 200 1200 B
python3 fw_rows.py refb EE 200
echo FW_QUEUE_B_DONE
