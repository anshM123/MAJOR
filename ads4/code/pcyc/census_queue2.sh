#!/bin/bash
# chain 2: exact-skeleton end-to-end check (pcyc_exact) then the extended census (pcyc_census).  Waits for the Delta = 100 D' chain.
ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$(dirname "$0")"
while pgrep -f "dq_prod.sh 100" > /dev/null; do sleep 20; done
python3 pcyc_exact.py 30 100 > cache/log_exact.txt 2>&1
for job in "300 4 600 17" "30 4 600 11" "100 4 600 11" "1000 4 600 11" "300 2 600 11" "300 8 600 11"; do
  set -- $job
  [ -f cache/census_D$1_nu$2_s$4.txt ] || python3 pcyc_census.py $1 $2 $3 $4 >> cache/census_queue.log 2>&1
done
