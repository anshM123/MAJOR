#!/bin/bash
# nf_pause.sh (NEARGAP_FILL.md): SIGSTOP lane A's python worker, run one short analysis command, SIGCONT it (keeps <= 2 processes).
# Usage: bash nf_pause.sh "cmd" logfile
P=$(pgrep -f "nf_lane.sh A" | head -1); W=$(pgrep -P "$P" -f python3 | head -1)
[ -n "$W" ] && kill -STOP "$W"
ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
cd "$(dirname "$0")"; eval "$1" > "$2" 2>&1; rc=$?
[ -n "$W" ] && kill -CONT "$W"
exit $rc
