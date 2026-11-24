#!/bin/sh
# graviton amplitude lanes (final code): $1 = C3H or EE; $2 = optional PID to wait for first
[ -n "$2" ] && while kill -0 $2 2>/dev/null; do sleep 15; done
ulimit -v 3000000
export OMP_NUM_THREADS=1
cd "$(dirname "$0")"
K=$1
S='0.015,0.02,0.025,0.03'
[ "$K" = C3H ] && S='0.015,0.02,0.025,0.03,0.035'
GC_DEG=6 python3 gm_spin.py $K 2 0.5,0.7,0.3 $S >> cache/run_spin_$K.log 2>&1
GC_DEG=6 python3 gm_spin.py $K 0 0.5,0.3,0.7 $S >> cache/run_spin_$K.log 2>&1
GC_DEG=6 python3 gm_spin.py $K 2 0.85 $S >> cache/run_spin_$K.log 2>&1
GC_DEG=6 python3 gm_spin.py $K 0 0.85 $S >> cache/run_spin_$K.log 2>&1
GC_DEG=6 python3 gm_spin.py $K 2 0.1 0.01,0.015 >> cache/run_spin_$K.log 2>&1
GC_DEG=6 python3 gm_spin.py $K 0 0.1 0.01,0.015 >> cache/run_spin_$K.log 2>&1
