#!/bin/bash
# HLJ_CLOSE: (364, 64) depth-2 engine in the fw (patched) basis (reference (2e4, 64) d2_p exists), then its normalisation.
cd "$(dirname "$0")"
ulimit -v 3000000
export BJ_ADAPT=1 BJ_PRUNE=1e-12 BJ_KTOL=1e-6 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 BJ_PATCH=1
python3 lj_rs_book.py 28 336 3 4 2 > cache/lj_rs_book_J28.log 2>&1
python3 lj_rs_eeq.py 364 64 2 4,5 fw > cache/lj_rs_eeq_D364_J64_d2.log 2>&1
