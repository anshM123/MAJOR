#!/bin/bash
# qb_lane.sh (QBIND_COROLLARY.md): run a list of qb_rows / qb_points jobs sequentially (one lane = one process at a time).
# Usage: bash qb_lane.sh <lane-name> "<args1>" "<args2>" ...   (each args string is passed to python3)
ulimit -v 3000000; export OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
cd "$(dirname "$0")"
lane=$1; shift
for job in "$@"; do
  echo "### $(date +%T) $job" >> cache/lane_$lane.log
  python3 $job >> cache/lane_$lane.log 2>&1
  echo "### $(date +%T) exit $?" >> cache/lane_$lane.log
done
