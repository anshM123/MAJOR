#!/bin/bash
# gs2_pause.sh CMD... : run a short analysis while one lane's compute process is SIGSTOPped (at most 2 processes compute at a time).
cd "$(dirname "$0")"
ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
P=${GS2_PID:-$(pgrep -f "^python3 gs2_rows.py" | sort -n | tail -1)}
[ -n "$P" ] && kill -STOP "$P"
trap '[ -n "$P" ] && kill -CONT "$P"' EXIT
"$@"
