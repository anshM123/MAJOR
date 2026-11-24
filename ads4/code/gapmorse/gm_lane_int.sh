#!/bin/sh
# xi-interval certificates (gm_cert8, 8 Taylor-model variables); run after lane $1 (a PID) has finished
[ -n "$1" ] && while kill -0 $1 2>/dev/null; do sleep 15; done
ulimit -v 3000000
export OMP_NUM_THREADS=1
cd "$(dirname "$0")"
for J in 0 7; do GC_DEG=8 python3 gm_cert8.py $J 0.5 0.005 0.02,0.025,0.03 >> cache/run_int.log 2>&1; done
GC_DEG=8 python3 gm_cert8.py 0 0.5 0.0025 0.025 >> cache/run_int.log 2>&1
GC_DEG=8 python3 gm_cert8.py 0 0.5 0.0125 0.02,0.025 >> cache/run_int.log 2>&1
for xi in 0.3 0.7; do GC_DEG=8 python3 gm_cert8.py 0 $xi 0.005 0.02,0.025 >> cache/run_int.log 2>&1; done
