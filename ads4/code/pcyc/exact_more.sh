#!/bin/bash
# more end-to-end scalar checks (D' vs exact skeleton) at Delta = 20, 50; runs after census_big.sh
ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$(dirname "$0")"
while pgrep -f "census_big.sh" > /dev/null; do sleep 10; done
for D in 20 50; do
  for g in 12,12,6,12,12,12,12 14,16,8,14,16,16,16; do
    [ -f cache/dq_D${D}_S_${g//,/-}.npz ] || python3 pcyc_dquad.py $D S $g > cache/log_dq_${D}_S_$g.txt 2>&1
  done
  python3 pcyc_exact.py $D 100 > cache/log_exact_$D.txt 2>&1
done
