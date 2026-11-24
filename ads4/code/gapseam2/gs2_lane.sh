#!/bin/bash
# gs2_lane.sh JOBFILE LOG : run the commands of JOBFILE one after the other (one compute process per lane; at most 2 lanes),
# ulimit -v 3 GB, one thread, no bytecode.  Each command's output -> its own log (redirect inside the job line).
cd "$(dirname "$0")"
ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
while IFS= read -r cmd; do
  [[ -z "$cmd" || "$cmd" == \#* ]] && continue
  echo "$(date +%H:%M:%S) start: $cmd" >> "$2"
  bash -c "$cmd"; rc=$?
  echo "$(date +%H:%M:%S) end ($rc): $cmd" >> "$2"
done < "$1"
