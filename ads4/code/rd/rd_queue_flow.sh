#!/bin/bash
# rd_flow production: two chains (<= 2 processes), each job ulimit -v 3 GB, 1 thread.
cd "$(dirname "$0")"
ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 RD_NSTEP=4000
chain() { for D in $1; do for nu in 2 4 6 8 12 16; do python3 rd_flow.py $D $nu 2048 2.0 > cache/rd_flow_D${D}_nu${nu}.log 2>&1; done; done; }
chain "30 1000" &
chain "100 300" &
wait
