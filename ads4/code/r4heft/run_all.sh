#!/bin/bash
# reproduces every number of R4ARC_HEFT.md (<= 2 processes at a time, 1 thread each, < 1 GB each)
cd "$(dirname "$0")"
ulimit -v 3000000; export OMP_NUM_THREADS=1
python3 rh_explore.py 2 > cache/explore_k2.txt 2>&1 & python3 rh_explore.py 4 > cache/explore_k4.txt 2>&1; wait
python3 rh_degree.py 2 > cache/degree_k2.log 2>&1 & python3 rh_degree.py 4 > cache/degree_k4.log 2>&1; wait
python3 rh_lemma_count.py > cache/lemmaY.log 2>&1
for k in 2 4; do
  RH_MODE=deg python3 rh_t2_heft.py run $k ab P1 > cache/t2d_P1_k${k}_ab.log 2>&1 &
  RH_MODE=deg python3 rh_t2_heft.py run $k ba P1 > cache/t2d_P1_k${k}_ba.log 2>&1; wait
  RH_MODE=deg python3 rh_t2_report.py $k P1 > /dev/null
done
python3 rh_t1_collinear.py scan 2 > cache/t1_scan_k2.log 2>&1 & python3 rh_t1_collinear.py scan 4 > cache/t1_scan_k4.log 2>&1; wait
python3 rh_t1_collinear.py shift 2 0.3 > cache/t1_shift_k2.log 2>&1 & python3 rh_t1_collinear.py shift 4 0.3 > cache/t1_shift_k4.log 2>&1; wait
python3 rh_t1_collinear.py report > /dev/null
python3 rh_t1_laurent.py 4 0.3 s D8 > cache/t1_laurent_k4_s.txt 2>&1 & python3 rh_t1_laurent.py 2 0.3 s D8 > cache/t1_laurent_k2_s.txt 2>&1; wait
python3 rh_t1_laurent.py 2 0.3 1243 D4 > cache/t1_laurent_k2_1243.txt 2>&1 & python3 rh_t1_tail.py 0.3 > /dev/null 2>&1; wait
python3 rh_t1_frozen.py 2 0.0375 > /dev/null 2>&1 & python3 rh_t1_frozen.py 4 0.0375 > /dev/null 2>&1; wait
RH_MODE=j python3 rh_t2_heft.py run 2 ab P1 > cache/t2_P1_k2_ab.log 2>&1 & RH_MODE=j python3 rh_t2_heft.py run 2 ba P1 > cache/t2_P1_k2_ba.log 2>&1; wait
RH_MODE=j python3 rh_t2_report.py 2 P1 > /dev/null
