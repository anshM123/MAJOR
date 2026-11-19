#!/bin/bash
# r1 queue A: (672, 36) K2 full nu-grid of the (336, 18) row (C2 + C3H, state + ref), then the reference basis for J = 36
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 r1_run.py k2 36 672 25.049,33.399,41.749,50.098,58.448,66.798,16.699,83.497,100.196,116.896,133.595,150.295,166.994,183.693,200.393,217.092,233.792,250.491,267.190
python3 r1_run.py refbasis 36
