#!/bin/bash
# pwimg (PW_IMAGES.md): reproduce everything, one process at a time (~1 CPU-h).
set -e
cd "$(dirname "$0")"
ulimit -v 3000000; export OMP_NUM_THREADS=1
PW_NSYM=8 python3 pw1_magic.py      # magic-block identity for the b-shifted blocks; exact (n <= 8) + 50-digit (n <= 40); e2e   ~8 min
python3 pw1b_e2e.py                 # refined end-to-end check against eestep/ee_quad direct quadrature                       ~10 min
python3 pw2_radial.py               # EE/scalar ratio c_s under radial-weight shifts (zeros at +-i/2 not protected)            ~0.5 min
python3 pw3_cost.py                 # exact cost of scalar-de-shift atoms, the monopole constant D                              ~15 min
python3 pw4_theoremW.py             # Theorem W without P3: R1' costs, K_tol, J*
python3 pw5_theoremL.py             # Theorem L: monopole tolerances and exact per-structure costs, J = 857 ... 1e6
python3 pw6_density.py              # closed-form densities of the index-shifted harmonics (Abel inversion), check
