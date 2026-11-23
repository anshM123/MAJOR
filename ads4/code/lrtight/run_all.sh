#!/bin/bash
# lrtight (LRINF_TIGHT.md): certified Lemma LR-inf with the tight constant int|psi|.  ~10 CPU-min, one process at a time.
set -e
cd "$(dirname "$0")"
ulimit -v 3000000; export OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
python3 lrt1_constants.py          # certified int|psi_k|, exact a1, a2           -> cache/lrt1_constants.*
python3 lrt2_m1.py                 # certified M1 lower bounds on beta-cells      -> cache/lrt2_m1.*   (~2 min)
python3 lrt3_rcert.py              # certified r(beta), rtilde, F_cert, tail      -> cache/lrt3_rcert.*
python3 lrt4_jstar.py              # J*(m) tables (Arb comparisons)              -> cache/lrt4_jstar.*
python3 lrt5_tolerances.py         # [DV] Theorem L' tolerances, sanity, hybrid  -> cache/lrt5_tolerances.*  (~4 min)
