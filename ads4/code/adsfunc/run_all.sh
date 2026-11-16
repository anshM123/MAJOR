#!/bin/bash
# Reproduces every number of ads4/ADS4_FUNCTIONAL.md (floating point; nothing here is a certificate).
# Outputs: cache/*.log (git-ignored).  Runtime: t1 ~1 min, t2 ~15-30 min, t3 ~1-2 h (machine-dependent), t4 ~1 min, t5 seconds.
set -e
cd "$(dirname "$0")"
mkdir -p cache
python3 -u t1_validate.py      | tee cache/t1.log
python3 -u t4_marginal.py      | tee cache/t4.log
python3 -u t5_light.py         | tee cache/t5.log
python3 -u t6_regge_estimate.py | tee cache/t6.log
python3 -u t2_regge_leading.py | tee cache/t2.log
python3 -u t3_hybrid_scan.py 20 50 100 300 1000 | tee cache/t3.log
