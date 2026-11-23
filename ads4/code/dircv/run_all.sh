#!/bin/bash
# DIRCV.md: reproduce the analysis (everything writes only to dircv/cache; dc_guard refuses other writes).
#   heavy (MC, ~2 CPU-h, two processes): ./run_new.sh (only its _x1/_x2 parts were used; _x3/_x4 were stopped), ./run_tail.sh
#   light (~5 CPU-min): below
set -e
cd "$(dirname "$0")"
ulimit -v 3000000; export OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
python3 dc_runall.py                     # pipelines for dir/cv x _w/_wP/_none; checks reproduction of pw_gamma.pkl
python3 dc_raw.py                        # raw state-level comparison (1100,100)
python3 dc_cvmean.py                     # the CV mean is exact (interpolation 1e-10, rho_2 total, reference rho_2)
python3 dc_compare.py                    # r with jackknife errors; paired same-CRN differences
python3 dc_sens.py                       # linear response of r to each node: dir - cv is the nu = 13.3 node
python3 dc_k3.py                         # k = 3 check with jackknife errors
python3 dc_scalseeds.py                  # scalar wave bias, all seeds pooled
python3 dc_mstest.py _x1,_x2             # lambda (graviton/scalar wave-bias ratio), 4 states + new seeds
python3 dc_pool.py 100 1100 _w,_x1,_x2   # pooled record set (8 super-replicas)
python3 dc_lambda.py _w
DC_LOCAL= python3 dc_lev.py _w
DC_LOCAL=_pool,_poolP python3 dc_lev.py _pool
python3 dc_tailan.py                     # per-point heavy-tail / CV diagnostics at nu = 13.3
DC_LOCAL= python3 dc_final.py _w
DC_LOCAL=_pool,_poolP python3 dc_final.py _pool
