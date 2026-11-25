#!/bin/bash
# N3C3H_FIX.md §3-4: validation (n1 vs patched n3, J = 26 and 29, x <= 0.93) and J = 26 edge nodes.  Two processes at most.
cd "$(dirname "$0")"
ulimit -v 3000000; export OMP_NUM_THREADS=1
printf '*\n!.gitignore\n' > cache/.gitignore
(python3 nf_run.py val_n1 n1 26 1026 0.85,0.87,0.90,0.92,0.93 > cache/val_n1_J26.out 2>&1
 python3 nf_run.py val_n3 n3 26 1026 0.85,0.87,0.90,0.92,0.93 > cache/val_n3_J26.out 2>&1
 python3 nf_run.py edge n3 26 1026 nu=966.664857,974.623870,981.866866,993.069476 6000 > cache/edge_n3_J26.out 2>&1) &
(python3 nf_run.py val_n1 n1 29 1029 0.85,0.88,0.90,0.92,0.93 > cache/val_n1_J29.out 2>&1
 python3 nf_run.py val_n3 n3 29 1029 0.85,0.88,0.90,0.92,0.93 > cache/val_n3_J29.out 2>&1) &
wait
