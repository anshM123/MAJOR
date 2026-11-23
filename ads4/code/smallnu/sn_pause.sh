#!/bin/bash
# sn_pause.sh: run a short analysis command while this note's sn_points workers are SIGSTOPped (keeps <= 2 computing processes).
# Usage: bash sn_pause.sh <command...>
P=$(pgrep -f "python3 sn_points.py" | head -1)
[ -n "$P" ] && kill -STOP $P
ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
cd "$(dirname "$0")"
"$@"; rc=$?
[ -n "$P" ] && kill -CONT $P
exit $rc
