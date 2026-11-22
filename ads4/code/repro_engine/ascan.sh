#!/bin/bash
# angular (phi, alpha) refinement scan at converged radial Y; state and reference in parallel (2 processes)
ulimit -v 3000000; export OMP_NUM_THREADS=1
cd "$(dirname "$0")"
for ns in 8,8,16,16,6,12,12 8,8,16,16,6,16,16 8,8,16,16,6,20,20 10,10,16,16,8,16,16; do
  python3 direct_rho.py 300 $ns f64 4.1 | tail -1 &
  python3 direct_rho.py 20000 $ns f64 4.1 | tail -1
  wait
done
echo ASCAN_DONE
