#!/bin/bash
# H_OO (HYPOTHESES_CLOSE2.md): K2 OO rows (3x3 p3 action, EE and C3H) on the lowest Delta_g = 300 rows of even J = 12..22 (the rows
# that carry the near-gap states), EE-row x grid; one cache file per state row.  Resumable.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for spec in "20 400" "22 336" "18 336" "18 400" "16 400"; do
  python3 oo_rows.py EE $spec; python3 oo_rows.py C3H $spec
done
export OO_XS=0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8
for spec in "14 400" "12 400"; do
  python3 oo_rows.py EE $spec; python3 oo_rows.py C3H $spec
done
echo "QUEUE OO-K DONE"
