#!/usr/bin/env bash
# Reproduces ads4/UNIFORM_REMAINDER.md.  Outputs in cache/ (git-ignored).  Uses <= 2 cores.
#   bash run_all.sh         everything (about 1.5 h on 2 cores; the e2e block code must be built: ads4/code/e2e/libcblock_*.so)
#   bash run_all.sh cert    the certificate only (a few minutes)
set -eu
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
mkdir -p cache
ANC=cache/regge_moments_order8.m
if [ ! -f "$ANC" ]; then
  echo "need the CHMRS ancillary file anc/regge_moments_order8.m of arXiv:2106.10274 at $ANC" >&2; exit 1
fi
python3 cert_regionA.py
[ "${1:-}" = "cert" ] && exit 0
python3 chmrs_series.py                      # parse the ancillary series (cached)
python3 v1_series_vs_exact.py > cache/v1.log # series vs exact blocks
python3 t1_structure.py                      # bulk-point structure, remainder coefficients
python3 residual_coeffs.py 3; python3 residual_coeffs.py 2   # curvature-subtracted remainder coefficients (both skeletons)
python3 t4_radius.py                         # Domb-Sykes estimate
python3 t2_bulkpoint_exact.py > cache/t2.log # exact blocks vs Regge / bulk-point at nu/m <= 0.5
python3 v2_mellin.py                         # independent J = 0 check (Mellin) vs series
UNIFORM_DF=3 python3 t3_scan.py 20 50 100 300 1000 > cache/t3_df3.log &
UNIFORM_DF=2 python3 t3_scan.py 20 50 100 300 1000 > cache/t3_df2.log
wait
python3 summarize.py > /dev/null
