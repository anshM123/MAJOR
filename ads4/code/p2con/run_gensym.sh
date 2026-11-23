#!/bin/bash
# generate the exact symbols at every Theorem-W sample state (2 processes)
cd "$(dirname "$0")"
ulimit -v 3000000; export OMP_NUM_THREADS=1
python3 pc_gensym.py $(cat cache/bpart_aa) > cache/gensym_a.log 2>&1 &
python3 pc_gensym.py $(cat cache/bpart_ab) > cache/gensym_b.log 2>&1 &
wait
