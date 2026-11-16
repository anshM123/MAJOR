#!/bin/sh
# Regenerates every number quoted in ads4/SPINNING_LOWK.md.   Usage:  sh run_all.sh
# Needs the e2e block-tensor cache ads4/code/e2e/data_modes/ (git-ignored; regenerate with  cd ../e2e && sh run_all.sh fast,  ~3.5 h).
# Caches of this directory: cache/ (git-ignored): lcmodes_*.npz (corner-completed mode data), data_t1_l*.pkl.
# Run time on 4 shared cores: corner completion ~10 min (l = 1) / ~20 min (l = 2) per case, v1 ~ 2 h, q1 ~ 1 h, c1 ~ 20 min.
set -e
cd "$(dirname "$0")"
gcc -O3 -march=native -fopenmp -shared -fPIC lcconv.c -o liblcconv.so
(cd ../e2e && gcc -O3 -march=native -fopenmp -shared -fPIC cblock.c -o libcblock_d.so -lm)
python3 c1_block_checks.py              # LC/tail block representation: H_0 vs 2F1, vs large-N radial, Casimir residual -> out_c1_block_checks.txt
python3 u1_uniformity.py                # uniformity of the heavy expansion of the shadow integrand  -> out_u1_uniformity.txt
python3 v1_scalar.py                    # scalar validation (corner-heavy scalar, k = 2, 3)        -> out_v1_scalar.txt
for b in 0.4 1.0 2.0; do python3 t1_lowk_e2e.py 1 build $b & done; wait   # corner completion (cached in cache/), 3 processes
for b in 0.4 1.0 2.0; do python3 t1_lowk_e2e.py 2 build $b & done; wait
python3 t1_lowk_e2e.py 1                # currents: completed moments vs leading formula            -> out_t1_lowk_e2e_l1.txt
python3 t1_lowk_e2e.py 2                # stress tensors                                            -> out_t1_lowk_e2e_l2.txt
python3 q1_grid.py                      # quadrature-grid check                                     -> out_q1_grid.txt
python3 v2_l2_validation.py             # LC vs radial at all converged points u < 0.25 (l = 1, 2)   -> out_v2_l2_validation.txt
python3 t2_corrections.py               # m^2(R-1), c1, c2, non-analytic-term fits                 -> out_t2_corrections.txt
python3 t3_tables.py                    # tables of SPINNING_LOWK.md                                -> out_t3_tables.txt
