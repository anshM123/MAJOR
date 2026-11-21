#!/bin/bash
# po_item3 (PARITY_ODD_CLOSE.md sect. 3): engine cross-check (po_k2cmp) then direct small-nu 3x3 nodes at even J = 6 (oo_points, unchanged)
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ulimit -v 3000000
python3 po_k2cmp.py > cache/po_k2cmp.log 2>&1
python3 oo_points.py 6 400 0.021,0.025,0.0293,0.0335,0.042,0.05 > cache/po_item3_oo_J6.log 2>&1
