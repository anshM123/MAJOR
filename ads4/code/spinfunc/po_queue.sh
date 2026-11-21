#!/bin/bash
# po_queue (PARITY_ODD_CLOSE.md sect. 2): parity-odd +- J = 5 EE k = 4 images on the analytic-radial engine (2 workers, <= 3 GB each).
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ulimit -v 3000000
case "$1" in
  A) python3 po_ee_odd.py run 20000 5 2 3 > cache/po_eeq_D20000_J5_K2.log 2>&1
     python3 po_ee_odd.py run 1005.0089775159122 5 4 5 > cache/po_eeq_D1005_J5_K4.log 2>&1 ;;
  B) python3 po_ee_odd.py run 305.0297506031432 5 4 5 > cache/po_eeq_D305_J5_K4.log 2>&1 ;;
esac
# second stage: the +- J = 5 binding-region state at Dg = 300 (m = 1.1935 Lam), for the m-dependence of s_d
if [ "$1" == "C" ]; then python3 po_ee_odd.py run 357.8 5 4 5 > cache/po_eeq_D357.8_J5_K4.log 2>&1; fi
