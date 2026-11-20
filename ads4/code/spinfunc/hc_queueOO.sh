#!/bin/bash
# H_OO: K2 rows of the even-J parity-odd (OO) block, EE and C3H, J = 4, 6, 8, 10, rows 400 / 1600 (resumable; re-run to resume)
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
for J in 4 6 8 10; do
  python3 hc_oo_rows.py EE $J
  python3 hc_oo_rows.py C3H $J
done
echo "QUEUE OO DONE"
