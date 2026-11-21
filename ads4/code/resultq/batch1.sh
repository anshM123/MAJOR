#!/bin/bash
# Result Q stage 0 LP scans (sequential; <= 2 threads)
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
cd "$(dirname "$0")"
for d in 4 6 10 12; do for dl in 0 0.05; do
  python3 lp_g4.py $d $d 3 --delta $dl | cut -c1-140
  python3 lp_g4.py $d $d 3 --delta $dl --matter | cut -c1-140
  python3 lp_g4.py $d 0 3 --delta $dl --no-b3 | cut -c1-140
done; done
for r in 2 5; do for dl in 0 0.05; do
  python3 lp_g4.py 8 8 $r --delta $dl | cut -c1-140
  python3 lp_g4.py 8 8 $r --delta $dl --psi3zero | cut -c1-140
  python3 lp_g4.py 8 0 $r --delta $dl --no-b3 | cut -c1-140
done; done
