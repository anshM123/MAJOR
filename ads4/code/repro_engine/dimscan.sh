#!/bin/bash
# per-dimension refinement scan (state and reference on the same grid)
ulimit -v 3000000; export OMP_NUM_THREADS=1
cd "$(dirname "$0")"
for ns in 12,8,8,8,6,8,8 8,12,8,8,6,8,8 8,8,12,8,6,8,8 8,8,8,12,6,8,8 8,8,8,8,9,8,8 8,8,8,8,6,12,8 8,8,8,8,6,8,12; do
  for D in 300 20000; do python3 direct_rho.py $D $ns f64 4.1 | tail -1; done
done
