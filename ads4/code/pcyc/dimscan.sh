#!/bin/bash
# per-dimension convergence scan of the D' rule (scalar): base grid + 6 nodes in one dimension at a time
ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$(dirname "$0")"
D=${1:-30}; K=${2:-S}
for g in 12,10,6,12,10,10,10 18,10,6,12,10,10,10 12,16,6,12,10,10,10 12,10,10,12,10,10,10 12,10,6,18,10,10,10 12,10,6,12,16,10,10 12,10,6,12,10,16,10 12,10,6,12,10,10,16; do
  [ -f cache/dq_D${D}_${K}_${g//,/-}.npz ] || python3 pcyc_dquad.py $D $K $g > cache/log_dq_${D}_${K}_$g.txt 2>&1
done
