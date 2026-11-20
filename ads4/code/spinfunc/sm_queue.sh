#!/bin/bash
# SMOOTH_BOUND.md §2: exact scalar remainder at near-gap states (sm_exact.py), sequential, 1 core.  Resumable (cached keys skipped).
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
SM_XPROBE=0.1:0.7125:0.025 SM_KEYTAG=d48 nice -n 5 python3 sm_exact.py 1000 48
nice -n 5 python3 sm_exact.py 1000 64
nice -n 5 python3 sm_exact.py 300 24,40
nice -n 5 python3 sm_exact.py 300 64
nice -n 5 python3 sm_exact.py 1000 100
nice -n 5 python3 sm_exact.py 300 100
echo "# queue done"
