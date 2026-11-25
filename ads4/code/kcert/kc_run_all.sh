#!/bin/bash
# K certificate (ads4/K_CERT.md): all outputs to cache/ ; ~1 CPU-min total, one process, one thread.
ulimit -v 3000000; export OMP_NUM_THREADS=1
cd "$(dirname "$0")"
python3 kc_dk0.py 300  > cache/kc_dk0_300.log
python3 kc_dk0.py 1000 > cache/kc_dk0_1000.log
python3 kc_grid.py     > cache/kc_grid.log
python3 kc_bounds.py   > cache/kc_bounds.log
python3 kc_selftest.py > cache/kc_selftest.log
