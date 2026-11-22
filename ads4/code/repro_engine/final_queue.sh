#!/bin/bash
# production + convergence/precision checks; at most 2 processes at a time
ulimit -v 3000000; export OMP_NUM_THREADS=1
cd "$(dirname "$0")"
A=10,10,16,16,8,16,16      # production grid
B=8,8,18,18,6,20,20        # independent second grid (convergence check to nu = 12)
C=8,8,16,16,6,16,16        # extended-precision check grid
python3 direct_rho.py 300 $A f64 12 | tail -1 &
python3 direct_rho.py 20000 $A f64 12 | tail -1; wait
python3 direct_rho.py 1000 $A f64 12 | tail -1 &
python3 direct_rho.py 300 $B f64 12 | tail -1; wait
python3 direct_rho.py 20000 $B f64 12 | tail -1 &
python3 direct_rho.py 1000 $B f64 12 | tail -1; wait
python3 direct_rho.py 300 $C f80 4.1 | tail -1 &
python3 direct_rho.py 20000 $C f80 4.1 | tail -1; wait
echo QUEUE_DONE
