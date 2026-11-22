#!/bin/bash
# cl_pause (RESULT_D_CLEANUP.md): run one analysis command while the edge-node worker (cl_queueA python if running, else cl_queueB/C) is SIGSTOPped, so that at most
# 2 of this agent's processes are busy; analyses are serialised by a lock.  usage: ./cl_pause.sh <command...>
cd "$(dirname "$0")"
exec 9> cache/cl_pause.lock
flock 9
P=$(pgrep -f "python3 cl_points.py 48 1048" | head -1)                                  # the J = 48 check (cl_queueC) first
[ -z "$P" ] && P=$(pgrep -f "python3 cl_points.py [0-9]+ [0-9]+ 0[.]0" | head -1)       # then the small-nu worker (cl_queueA)
[ -z "$P" ] && P=$(pgrep -f "python3 cl_points.py [0-9]+ [0-9]+ 0[.][1-9]" | head -1)   # else the edge / J = 48 worker
[ -n "$P" ] && kill -STOP $P
ulimit -v 3000000
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$@"; rc=$?
[ -n "$P" ] && kill -CONT $P
exit $rc
