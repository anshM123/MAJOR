#!/bin/bash
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
while kill -0 17915 2>/dev/null; do sleep 20; done
python3 y1_rows.py evenx 8,6 1.1,1.2,1.3
