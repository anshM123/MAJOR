#!/bin/bash
# N3C3H_FIX.md §4: the rest of the J = 26 edge (GK nodes 0.9357, 0.9746), the n1 basis anchors at two batch-3 n1 nodes, and n3 at the
# two batch-3 overlap nodes (0.9198, 0.9276).  Started after run_val.sh's edge jobs; at most two processes.
cd "$(dirname "$0")"
ulimit -v 3000000; export OMP_NUM_THREADS=1
while pgrep -f "nf_run.py edge2 " > /dev/null; do sleep 20; done
python3 nf_run.py edge3 n3 26 1026 nu=958.315057,998.156510 6000 > cache/edge3_n3_J26.out 2>&1 &
while pgrep -f "nf_run.py edge n3" > /dev/null; do sleep 20; done
python3 nf_run.py anchor n1 26 1026 nu=942.006243,877.94542 > cache/anchor_n1_J26.out 2>&1
python3 nf_run.py ovl n3 26 1026 nu=942.006243,949.965255 6000 > cache/ovl_n3_J26.out 2>&1
wait
