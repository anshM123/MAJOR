#!/bin/bash
# Route W (LARGESPIN_W.md): reproduce everything (~0.6 CPU-h, one process).
set -e
ulimit -v 3000000; export OMP_NUM_THREADS=1
cd "$(dirname "$0")"
python3 lsw_weyl.py        # §2.2  Weyl-reflection parity check          -> cache/weyl.txt
python3 lsw_costs.py       # §4.1  per-unit Lemma-I costs (slow part)     -> cache/costs_1000.{json,txt}
python3 lsw_theorem.py     # §4.1  sign-free uniform-K tolerances          -> cache/theorem_1000.{json,txt}
python3 lsw_classwise.py   # §4.2  class-wise tolerances, scalar multiples -> cache/classwise_1000.txt
python3 lsw_jstar.py       # §4.3  J* for stated premise constants         -> cache/jstar_1000.txt
python3 lsw_r1rank.py      # §9.1  rank test of the leading-r EE kernel   -> cache/r1rank.txt
