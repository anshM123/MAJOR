#!/bin/bash
# H_OO (HYPOTHESES_CLOSE2.md): direct small-nu n3 nodes of the 3x3 [even, odd] action at even J (OO block), Delta_g = 300 rows
# carrying the near-gap states; state + reference (2e4) on CRN, one cache file per state.  Resumable (finished nodes skipped).
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OO_M=10 OO_R=8
X336=0.025,0.03,0.035,0.04,0.05,0.0625,0.075,0.1
X400=0.021,0.025,0.0293,0.0335,0.042,0.05,0.075,0.1
python3 oo_points.py 18 336 $X336
python3 oo_points.py 18 400 $X400
python3 oo_points.py 20 400 $X400
python3 oo_points.py 22 336 $X336
python3 oo_points.py 16 400 $X400
python3 oo_points.py 14 400 $X400
python3 oo_points.py 12 400 $X400
echo "QUEUE OO-D DONE"
