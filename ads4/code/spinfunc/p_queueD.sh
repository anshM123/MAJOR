#!/bin/bash
# PARITY_ODD_ROUTEA.md sect. 11: direct (RQMC) odd-J rows at Delta_g = 100, kind C2 (resumable: re-run the script).
cd /home/user/MAJOR/ads4/code/spinfunc
export OMP_NUM_THREADS=1 L1_DEG=6
while pgrep -f "p8_odd_pw.py validate" > /dev/null; do sleep 15; done
P8_ROWS=100 python3 p8_odd_pw.py run C2 5 7 9 >> cache/p8_run_C2.log 2>&1
P8_ROWS= P8_MFL=1.05 python3 p8_odd_pw.py run C2 5 >> cache/p8_run_C2.log 2>&1
P8_ROWS= P8_MFL=1.1,1.2,1.3 P8_MHI=13 python3 p8_odd_pw.py run C2 5 7 9 >> cache/p8_run_C2.log 2>&1
P8_MFL=1.05 python3 p8_k2rows.py C2 5 >> cache/p8_k2rows_C2.log 2>&1
P8_MFL=1.1,1.2,1.3 python3 p8_k2rows.py C2 5 7 9 >> cache/p8_k2rows_C2.log 2>&1
