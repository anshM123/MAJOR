#!/bin/bash
# ng_pause (NEARGAP2.md): run a short analysis while SIGSTOPping the lane-B worker (so that at most two of this note's processes are
# busy at any time; the pattern of spinfunc/cl_pause.sh).  usage: ./ng_pause.sh <command...>
P=$(pgrep -f "ng_rows.py|ng_points.py|ng_engine5.py" | head -1)
[ -n "$P" ] && kill -STOP $P
ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
"$@"; rc=$?
[ -n "$P" ] && kill -CONT $P
exit $rc
