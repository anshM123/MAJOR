#!/bin/sh
# GAP_C3 reproduction (~1.5 CPU-h; at most 2 processes).  Outputs in cache/.
ulimit -v 3000000
export OMP_NUM_THREADS=1
cd "$(dirname "$0")"
# 1. certified scalar nodes (log-factored integrand, Taylor models of degree 6), J x xi x rho
GC_DEG=6 python3 gc_run.py 0,2,4,5 0.5,0.3,0.7,0.1,0.85 0.01,0.015,0.02,0.025 > cache/run_a.log 2>&1 &
GC_DEG=6 python3 gc_run.py 6,7,8 0.5,0.3,0.7,0.1,0.85 0.01,0.015,0.02,0.025 > cache/run_b.log 2>&1
wait
# 2. degree-8 check at xi = 0.5
GC_DEG=8 python3 gc_run.py 0 0.5 0.02,0.025,0.03 > cache/run_c.log 2>&1 &
# 3. spinning amplitudes (project builders on the Taylor-model backend), single-structure J
GC_DEG=5 sh -c "python3 gc_spin.py EE 2 0.3 0.01,0.015,0.02; python3 gc_spin.py EE 2 0.5,0.7 0.01,0.015; \
  python3 gc_spin.py C3H 2 0.5,0.7 0.01,0.015; python3 gc_spin.py EE 0 0.5 0.01,0.015; python3 gc_spin.py C3H 0 0.5 0.01,0.015" > cache/run_spin.log 2>&1
wait
# 4. [NE] ideal (sampled) sup norms
python3 gc_ideal.py S 0 0.5 0.03,0.05,0.1 > cache/run_ideal.log 2>&1
python3 gc_ideal.py S 7 0.5 0.03,0.05,0.1 >> cache/run_ideal.log 2>&1
python3 gc_ideal.py EE 2 0.5 0.03,0.05,0.1 >> cache/run_ideal.log 2>&1
# 5. tables
python3 gc_table.py > /dev/null
python3 gc_report.py > cache/gc_report.txt
