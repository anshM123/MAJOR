#!/usr/bin/env bash
# Reproduces ads4/OMEGA_U.md.  Outputs in cache/ (git-ignored).  Uses <= 2 cores.
#   bash run_all.sh          everything (~1.5-2 h on 2 cores; the u4 scans dominate)
#   bash run_all.sh fast     u1, u2, u3, u5, u6 only (~5 min; u3 builds the CHMRS fixed-eta coefficient caches on first run)
# Needs (read-only): ../uniform/cache/chmrs_series.pkl and ../uniform/cache/remainder_coeffs.pkl (made by ../uniform/run_all.sh from the
# CHMRS ancillary file), ../robust/functionals_certified.json, ../robust/{rcommon,check,lp}.py.
set -eu
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
mkdir -p cache
python3 u1_flat.py   > /dev/null          # exact + numeric flat-limit identity, amplification, odd-prefactor sensitivity -> cache/u1.log
python3 u2_coeffs.py > /dev/null          # remainder coefficients k = 2..7, Df = 3, 2, 2.5                               -> cache/u2.log
OMEGAU_DF=3 python3 u3_flatlimit.py > /dev/null   # AdS skeleton action -> flat action as Dg -> infinity             -> cache/u3_Df3.log
python3 u5_light.py  > /dev/null          # light side, leverages, bound vs n                                            -> cache/u5.log
python3 u6_regge.py  > /dev/null          # leading Regge action, Omega^U vs Omega                                       -> cache/u6.log
[ "${1:-}" = "fast" ] && exit 0
OMEGAU_DF=3 python3 u4_scan.py 100 20 50 300 1000 > cache/u4_run_Df3.log &
OMEGAU_DF=2 python3 u4_scan.py 100 20 50 300 1000 > cache/u4_run_Df2.log
wait
python3 u7_summary.py > /dev/null          # compact tables -> cache/u7_summary.txt
