#!/bin/bash
# hs_run (HIGHSPIN.md): run one short analysis while the most recently started hs_points worker is SIGSTOPped, so that at most two
# processes of this note compute at any time.  Usage: ./hs_run.sh 'cmd'
cd "$(dirname "$0")"
p=$(pgrep -n -f "python3 hs_points.py")
[ -n "$p" ] && kill -STOP $p
( ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1; eval "$1" ); rc=$?
[ -n "$p" ] && kill -CONT $p
exit $rc
