#!/bin/bash
# r1 queue B: direct (RQMC, scalar-row-free point_g) C2 at (672, 36) on the six direct nodes of the (336, 18) row, state + ref
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 r1_run.py direct C2 36 672 33.399,25.049,41.749,50.098,58.448,66.798
