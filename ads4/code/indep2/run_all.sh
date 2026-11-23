#!/bin/bash
# Independent check 2 (INDEP_CHECK2.md). About 1.5 CPU-h in total; at most 2 processes. Outputs in cache/.
set -e
cd "$(dirname "$0")"
ulimit -v 3000000; export OMP_NUM_THREADS=1
python3 -B t1_magic.py            > cache/t1_magic.log            # ~30 min (exact sympy, 50-digit, controls, pointwise)
python3 -B t1b_monopole.py        > cache/t1b_monopole.log        # ~8 min
python3 -B t1c_even_analytic.py   > cache/t1c_even_analytic.log   # ~1 min
python3 -B t1d_pointwise.py       > cache/t1d_pointwise.log       # ~10 min
python3 -B t2_lrinf.py            > cache/t2_lrinf.log            # ~8 min (beta grid 1..32)
python3 -B t2_lrinf.py 0.3 0.5 0.7 0.85 > cache/t2_lrinf_extra.log
python3 -B t2b_jstar_m.py         > cache/t2b_jstar_m.log
python3 -B t2c_rbound.py          > cache/t2c_rbound.log          # ~20 min (stated bound on the stated grid)
python3 -B t2c_rbound.py 0.5 0.6 0.7 0.8 0.9 0.95 > cache/t2c_rbound_extra.log
python3 -B t4_theoremW.py         > cache/t4_theoremW.log         # ~2 min (reads code/lsw/cache/costs_1000.json as data)
python3 -B t5_counts.py           > cache/t5_counts.log
