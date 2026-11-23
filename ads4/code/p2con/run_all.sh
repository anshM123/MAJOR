#!/bin/bash
# P2_CONTRACT.md: reproduce everything (~4 CPU-h; at most 2 processes).  $1 = path to CHMRS ancillary regge_moments_order8.m.
set -e
cd "$(dirname "$0")"
ulimit -v 3000000; export OMP_NUM_THREADS=1
printf '*\n!.gitignore\n' > cache/.gitignore
python3 -c "import json; R=json.load(open('../lsw/cache/costs_1000.json')); open('cache/betas.txt','w').write('\n'.join(repr(r['beta']) for r in R)+'\n')"
split -n l/2 cache/betas.txt cache/bpart_
./run_gensym.sh                          # §1.2 exact symbols at all 47 states (~1.6 h wall, 2 proc)
python3 pc_checks_ksh.py                 # §1.7 O(r) structure checks
python3 pc_run.py 0 24 & python3 pc_run.py 24 47; wait          # §1.3-1.6 block side (tables cached first call)
python3 pc_ksh_run.py 0 24 & python3 pc_ksh_run.py 24 47; wait  # §1.7 kernel O(r) term
python3 pc_p4.py; python3 pc_flat.py "$1"; python3 pc_p4b.py "$1"   # §2
for m in full block bh; do PC_MODE=$m python3 pc_jstar.py; done   # §4
