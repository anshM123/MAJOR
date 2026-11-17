#!/bin/bash
# SUPERCONVERGENT_UPLIFT.md: exact (sympy) leading-Regge responses from ads4/code/e2e/laplace.py M-tables (imported read-only). Single core, ~10 min.
set -e
set -o pipefail
cd "$(dirname "$0")"
mkdir -p cache
nice python3 -u s0_check_D1sq.py | tee cache/s0.log
nice python3 -u s1_responses.py  | tee cache/s1.log
nice python3 -u s2_even.py       | tee cache/s2.log
nice python3 -u c1_flip_contamination.py | tee cache/c1.log
nice python3 -u c2_sc_flip.py    | tee cache/c2.log
