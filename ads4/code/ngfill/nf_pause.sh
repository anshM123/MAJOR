#!/bin/bash
# nf_pause.sh (NEARGAP_FILL.md): SIGSTOP one running compute worker (default: the newest 'python3 nf_points.py' process; override with
# PAUSE_PAT), run one short analysis command, SIGCONT it -- keeps this note at <= 2 computing processes.
# Usage: [PAUSE_PAT=...] bash nf_pause.sh "cmd" logfile
W=$(pgrep -n -f "${PAUSE_PAT:-python3 nf_points.py}")
[ -n "$W" ] && kill -STOP "$W"
ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
cd "$(dirname "$0")"; eval "$1" > "$2" 2>&1; rc=$?
[ -n "$W" ] && kill -CONT "$W"
exit $rc
