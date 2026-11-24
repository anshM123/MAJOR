#!/bin/bash
# GAP_UNIFORM.md reproduction (~0.6 CPU-h, <= 2 processes).  Reads gapnear/cache (read only); writes only gapuni/cache.
set -e
cd "$(dirname "$0")"
ulimit -v 3000000
export OMP_NUM_THREADS=1
python3 gu1_tol.py 60 > cache/gu1_run.log            # ~4 min: Lemma T check, tolerances J = 0..60 at 1000/2000/5000
python3 gu1b_weyl.py                                  # ~1 min: 2x2 Weyl tolerances (even J)
python3 gu2_growth.py                                 # seconds: Gamma rows 1000/1600/2000 (needs gapnear/cache/gn_hmono_*.json)
python3 gu2_growth.py inc                             # ~1 min: incremental criterion (appends to gu2_growth.txt)
python3 gu3_seam.py                                   # ~2 min: seam tolerances
python3 gu5_cover.py                                  # seconds: growth-allowance table
rm -f cache/gu4_scalar.txt
python3 gu4_scalar.py 5:1.05,12:1.05 1000,1500,2000,3000,5000,10000 > cache/gu4_runA.log 2>&1 &
python3 gu4_scalar.py 26:1.05,8:1.3 1000,1500,2000,3000,5000,10000 > cache/gu4_runB.log 2>&1   # ~15 min each
wait
