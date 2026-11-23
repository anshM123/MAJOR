#!/bin/bash
# sn_run.sh: one compute lane (sequential jobs from a job file), memory/thread limits of SMALLNU_WINDOW.md.
# Usage: bash sn_run.sh jobfile logfile      (each line of jobfile: arguments of sn_points.py)
ulimit -v 3000000
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
cd "$(dirname "$0")"
while read -r line; do
  [ -z "$line" ] && continue
  case "$line" in \#*) continue;; esac
  echo "## $(date +%T) start: $line" >> "$2"
  python3 sn_points.py $line >> "$2" 2>&1
  echo "## $(date +%T) end (rc=$?): $line" >> "$2"
done < "$1"
