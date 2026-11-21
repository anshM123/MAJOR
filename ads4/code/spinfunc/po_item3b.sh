#!/bin/bash
# po_item3b (PARITY_ODD_CLOSE.md sect. 3): direct small-nu nodes at odd J = 15 (hc_points n3 path, unchanged), Delta = 400 row
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ulimit -v 3000000
R=${1:-400}; XS=0.021,0.025,0.0293,0.0335,0.042,0.05; [ "$R" == "336" ] && XS=0.025,0.03,0.035,0.04,0.05,0.0625; Z2_NUSPLIT=5 HC_METHOD=n3 python3 hc_points.py 15 $R $XS > cache/po_item3_hc_J15_$R.log 2>&1
# second stage (row 336, where the near-gap odd J = 15 states are held): ./po_item3b.sh 336
