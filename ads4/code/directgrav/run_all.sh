#!/bin/bash
# Direct graviton near-gap test (DIRECT_GRAVITON_NEARGAP.md).  <= 2 cores.  Caches/logs in cache/ (git-ignored).
# Reads (read-only): fourpt/tttt (TaylorWaves, t13/t16 prediction machinery), e2e, exactskel (ExactState, e1/e3/e7 caches), omegau.
set -e
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
mkdir -p cache
python3 c1_calib.py 100 1.04 4 3                       # position-space quadrature vs exact skeleton, usable nu range (~5 min)
python3 gwaves.py 100 1.04 4 > cache/gw_A.log 2>&1 &   # exact <TTTT> waves, ~2.7 h each (resumable)
python3 gwaves.py 100 1.3 6  > cache/gw_C.log 2>&1
wait
python3 a1_moments.py cache/gw_100_1.040_4.npz 100 1.04 16 > cache/a1_A.log 2>&1 &
python3 a1_moments.py cache/gw_100_1.300_6.npz 100 1.3 16  > cache/a1_C.log 2>&1
wait
rm -f cache/a2_summary.txt
python3 a2_action.py cache/a1_gw_100_1.040_4.pkl cache/a1_gw_100_1.300_6.pkl
