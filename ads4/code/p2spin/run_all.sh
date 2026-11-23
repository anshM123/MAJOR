#!/bin/bash
# P2_SPINNING.md: reproduce everything (~0.5 CPU-h, one process).
set -e
ulimit -v 3000000; export OMP_NUM_THREADS=1
cd "$(dirname "$0")"
python3 s1_structures.py   # §1.1 conserved <TTO_{Delta,J}> over Q(Delta,J), counts       -> cache/s1_structures.{txt,pkl}
python3 s1b_canonical.py   # §1.2 canonical exact bases                                     -> cache/s1b_canonical.{txt,pkl}
python3 s1c_validate.py    # §1.3 integer J, J = 0 vs tto_conserved, mpmath conservation    -> cache/s1c_validate.txt
python3 s2_main.py         # §2   heavy symbols to O(1/m^2), helicity frame, reflection     -> cache/s2_main.{txt,pkl}
python3 s3_spin.py         # §3   CHMRS mechanism: seed validation + helicity-diagonal model -> cache/s3_spin.txt
python3 s3b_leakage.py     # §3.3 helicity leakage sizes                                     -> cache/s3b_leakage.txt
python3 s4_jstar.py        # §4   Theorem W J* with class-wise constants                     -> cache/s4_jstar.txt
