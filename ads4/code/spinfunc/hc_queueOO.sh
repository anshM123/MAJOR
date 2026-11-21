#!/bin/bash
# H_OO: K2 rows of the even-J parity-odd (OO) block, EE and C3H, rows 400 / 1600 (resumable; re-run to resume).
# J = 4, 6: 19 nodes (done first); J = 8, 10: 11 nodes (cost).
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
while pgrep -f "hc_oo_rows.py" > /dev/null; do sleep 30; done
for J in 4 6; do python3 hc_oo_rows.py EE $J; python3 hc_oo_rows.py C3H $J; done
export HC_XS=0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8
for J in 8 10; do python3 hc_oo_rows.py EE $J; python3 hc_oo_rows.py C3H $J; done
echo "QUEUE OO DONE"
