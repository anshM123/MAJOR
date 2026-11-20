#!/bin/bash
# z2 queue S2: independent repeat (new seed, 2x points) of the J = 28 thimble nodes nu = 5.0, 6.67 (tagged cache z2s2_*), reproducibility test
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 Z2_TAG=s2 Z2_SEED=5151 Z2_M=10
python3 z2_points.py 28 4.999254944,6.665673259 336
