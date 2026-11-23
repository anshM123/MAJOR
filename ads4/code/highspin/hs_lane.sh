#!/bin/bash
# hs_lane (HIGHSPIN.md): run a list of hs_* jobs sequentially in ONE process lane (at most two lanes at once), each under
# ulimit -v 3000000 and one thread; log to cache/lane_<name>.log.  Usage: ./hs_lane.sh name 'cmd1' 'cmd2' ...
cd "$(dirname "$0")"
name=$1; shift
for c in "$@"; do
  echo "$(date +%H:%M:%S) start: $c" >> cache/lane_$name.log
  ( ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1; eval "$c" ); rc=$?
  echo "$(date +%H:%M:%S) end ($rc): $c" >> cache/lane_$name.log
done
