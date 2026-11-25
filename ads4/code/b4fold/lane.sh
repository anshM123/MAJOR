#!/bin/bash
# b4fold compute/analysis lane: ulimit -v 3000000, one thread, no bytecode; each command logged with times.  Usage: bash lane.sh LOG cmd1 ';;' cmd2 ...
ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
cd "$(dirname "$0")"
LOG="cache/$1"; shift
cmd=""
run() { echo "### $(date -u +%H:%M:%S) $cmd" >> "$LOG"; eval "$cmd" >> "$LOG" 2>&1; echo "### $(date -u +%H:%M:%S) exit $?" >> "$LOG"; }
for a in "$@"; do if [ "$a" = ";;" ]; then run; cmd=""; else cmd="$cmd $a"; fi; done
[ -n "$cmd" ] && run
