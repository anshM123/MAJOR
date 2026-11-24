#!/bin/sh
# GAP_MORSE reproduction (~5.5 CPU-h; at most 2 processes; outputs in cache/).  Lanes A and B run in parallel.
ulimit -v 3000000
export OMP_NUM_THREADS=1
cd "$(dirname "$0")"
R='0.03,0.04,0.05,0.06,0.07'; R1='0.02,0.025,0.03,0.035,0.04'
# 1. scalar k = 2 integrand with Q_J: certified Taylor-model data on Hermitian balls, degree 8 (gm_cert via gm_run)
( GC_DEG=8 python3 gm_run.py 0,2,4,5 0.5,0.3,0.7,0.85 $R; GC_DEG=8 python3 gm_run.py 0,2,4,5 0.1 $R1 ) > cache/run_a.log 2>&1 &
( GC_DEG=8 python3 gm_run.py 6,7,8 0.5,0.3,0.7,0.85 $R;   GC_DEG=8 python3 gm_run.py 6,7,8 0.1 $R1 ) > cache/run_b.log 2>&1
wait
# 2. degree-10 phase-only records (J-independent; smaller remainder of the phase -> smaller kappa)
( GC_DEG=10 python3 gm_cert.py phase 0.1 0.03,0.035,0.04; GC_DEG=10 python3 gm_cert.py phase 0.3 0.05,0.06;
  GC_DEG=10 python3 gm_cert.py phase 0.5 0.06,0.07 ) > cache/run_phase.log 2>&1 &
( GC_DEG=10 python3 gm_cert.py phase 0.85 0.05,0.06; GC_DEG=10 python3 gm_cert.py phase 0.7 0.06,0.07 ) > cache/run_phase2.log 2>&1
wait
# 3. graviton amplitudes EE (k = 2) and C3H (k = 3), J = 0, 2 (degree 6; frozen <TTO> vectors at Delta = 1002)
sh gm_laneS.sh EE &
sh gm_laneS.sh C3H
wait
# 4. xi-interval certificates (8 Taylor-model variables)
sh gm_lane_int.sh
# 5. tables
python3 gm_table.py > cache/run_table.log 2>&1
python3 gm_table.py int > cache/run_int_table.log 2>&1
python3 gm_report.py > /dev/null
