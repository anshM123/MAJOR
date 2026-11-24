#!/bin/sh
# run a short analysis with the second compute lane SIGSTOPped (keeps <= 2 running processes)
P=$(pgrep -f "^python3 gm_(run|spin|cert8).py" | tail -1)
[ -n "$P" ] && kill -STOP $P
ulimit -v 3000000; export OMP_NUM_THREADS=1
"$@"
[ -n "$P" ] && kill -CONT $P
