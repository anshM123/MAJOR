#!/bin/bash
# larger census: 4000 starts per (Delta, nu); one chain; runs after census_queue2.sh
ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$(dirname "$0")"
while pgrep -f "census_queue2.sh" > /dev/null; do sleep 10; done
for D in 30 100 300 1000; do for nu in 2 4 8; do
  [ -f cache/census_D${D}_nu${nu}_s23.txt ] || python3 pcyc_census.py $D $nu 4000 23 >> cache/census_big.log 2>&1
done; done
