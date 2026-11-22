#!/bin/bash
# radial-Y refinement scan; state and reference in parallel (2 processes)
ulimit -v 3000000; export OMP_NUM_THREADS=1
cd "$(dirname "$0")"
for ny in 16 20 24; do
  ns=8,8,$ny,$ny,6,8,8
  python3 direct_rho.py 300 $ns f64 4.1 | tail -1 &
  python3 direct_rho.py 20000 $ns f64 4.1 | tail -1
  wait
done
