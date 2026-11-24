#!/bin/bash
# EDGE_R.md: reproduce (~1.9 CPU-h; at most 2 processes).  Read-only imports from highspin/, agproof/, p4win/, lsw/, spinfunc/, exactskel/.
set -e
cd "$(dirname "$0")"
ulimit -v 3000000; export OMP_NUM_THREADS=1
printf '*\n!.gitignore\n' > cache/.gitignore
python3 er_diag.py 29 1029 0.92,0.94,0.95,0.96,0.97 > cache/diag.out          # §2.1 saddle / masks
python3 er_diag2.py 29 1029 0.93 10 > cache/diag2_093.out; python3 er_diag2.py 29 1029 0.96 10 > cache/diag2_096.out   # §2.1 proposals
python3 er_cut.py 29 1029 0.90,0.93,0.95,0.96,0.97 12 1.0 0,3,4,5,6,8 > cache/cut.out  # §2.1 radial cuts (no plateau)
python3 er_n3test.py 29 1029 0.93,0.96 9 3.5 > cache/n3test.out                  # §2.2 n3 scalar
python3 er_scalar.py 26,30,34,38,40,48 > cache/scalar.out &                      # §4 exact scalar
python3 er_sval.py 29 0.88,0.90,0.92,0.935,0.95,0.96,0.967,0.972 > cache/sval_J29.out; wait   # §2.3 n3 vs exact scalar
ER_M=10 python3 er_points.py 29 1029 0.936549 > cache/pilot_n3_J29_a.out &       # §2.4 n3 graviton pilots
ER_M=10 python3 er_points.py 29 1029 0.9605443 > cache/pilot_n3_J29_b.out; wait
ER_METHOD=n1 ER_M=11 python3 er_points.py 26 1026 0.93 > cache/pilot_J26_n1.out
ER_SUB=n3ov ER_M=10 python3 er_points.py 26 1026 0.93 > cache/pilot_J26_n3a.out &
ER_M=10 python3 er_points.py 26 1026 0.96 > cache/pilot_J26_n3b.out; wait
python3 er_n3edge.py > cache/n3edge.out                                          # §2.4 Gamma at the n3 nodes
ER_NO_N3=1 python3 er_row.py 29 > /dev/null                                       # §3 J = 29 (EDGE-R) budget, n1 data (n3 C3H unusable)
python3 er_jobs.py                                                               # §6 ext/jobs_batch3.txt
