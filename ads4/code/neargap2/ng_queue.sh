#!/bin/bash
# ng_queue (NEARGAP2.md): a job lane.  Waits for pid $2 (if given), then pops and runs the first line of cache/queue_$1.txt until the
# file is empty; each line is a shell command run from neargap2/ under ulimit -v 3 GB, one thread.  Two lanes at most (A, B).
cd "$(dirname "$0")"
ulimit -v 3000000
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
[ -n "$2" ] && while kill -0 $2 2>/dev/null; do sleep 15; done
Q=cache/queue_$1.txt
while true; do
  [ -s $Q ] || { sleep 30; [ -s $Q ] || { [ -e cache/queue_$1.stop ] && exit 0; continue; }; }
  L=$(head -n 1 $Q); sed -i '1d' $Q
  echo "$(date +%H:%M:%S) start: $L" >> cache/lane_$1.log
  bash -c "$L"
  echo "$(date +%H:%M:%S) end ($?): $L" >> cache/lane_$1.log
done
