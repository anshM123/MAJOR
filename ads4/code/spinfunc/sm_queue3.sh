#!/bin/bash
# SMOOTH_BOUND.md §2.2: Delta_g = 300 near-gap exact scalar with N0 = 20-60 (sm_nsum: stable to 1e-10 at (400, 100), N0 >= 80 is not).
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
SM_NCONF=20:960:1,60:960:1,40:960:1 SM_KEYTAG=n0t nice -n 5 python3 sm_exact.py 300 24,40,64,100
SM_NCONF=20:960:1,60:960:1,40:960:1 SM_KEYTAG=n0t nice -n 5 python3 sm_exact.py 1000 100
echo "# queue3 done"
