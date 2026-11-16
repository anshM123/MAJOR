#!/bin/bash
# Reproduces every number of ads4/ADS4_LIGHT.md (floating point / sympy exact; nothing here is a certificate).
# Outputs: cache/*.log (git-ignored).  Runtime: t1 ~ 4-5 min (mpmath Mellin integrals), t2 ~ 1 min, t3, t4 seconds.
set -e
set -o pipefail
cd "$(dirname "$0")"
mkdir -p cache
python3 -u t2_flat_split.py    | tee cache/t2.log
python3 -u t3_light_omega.py   | tee cache/t3.log
python3 -u t4_sensitivity.py   | tee cache/t4.log
python3 -u t1_scalar_mellin.py | tee cache/t1.log
