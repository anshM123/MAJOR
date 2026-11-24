#!/bin/bash
# nf_lane.sh (NEARGAP_FILL.md): run a list of commands sequentially under the resource limits; log to cache/lane_<name>.log
# Usage: bash nf_lane.sh name "cmd1" "cmd2" ...
ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
cd "$(dirname "$0")"
name=$1; shift
for c in "$@"; do
  echo "### $(date +%T) $c" >> cache/lane_$name.log
  eval "$c" >> cache/lane_$name.log 2>&1
  echo "### $(date +%T) exit $?" >> cache/lane_$name.log
done
