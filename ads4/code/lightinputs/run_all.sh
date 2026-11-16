#!/bin/bash
# Reproduces the computations of ads4/LIGHT_INPUTS.md (sympy exact / numpy; single core; < 1 min total).
# Outputs: cache/*.log (git-ignored).
set -e
set -o pipefail
cd "$(dirname "$0")"
mkdir -p cache
python3 -u t1_frame_regularity.py  | tee cache/t1.log
python3 -u t2_shock_commutator.py  | tee cache/t2.log
python3 -u t3_clsd_counting.py     | tee cache/t3.log
python3 -u t4_mtable_poles.py      | tee cache/t4.log
