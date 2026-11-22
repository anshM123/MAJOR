#!/bin/bash
# HLJ_CLOSE §3.3: (336, 32) at depth 3 (relative order a^2), then its validation (out-of-sample test of the J = 28 tail factors) and window chain.
cd "$(dirname "$0")"
ulimit -v 3000000
export BJ_ADAPT=1 BJ_PRUNE=1e-12 BJ_KTOL=1e-6 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
while pgrep -f "lj_rs_queueR.sh" > /dev/null; do sleep 30; done
python3 lj_rs_eeq.py 336 32 3 4 win > cache/lj_rs_eeq_D336_J32_d3.log 2>&1
python3 lj_rs_val.py 32 3 4 2 > cache/lj_rs_val_J32_d3.log 2>&1
source ./lj_rs_env32.sh
L=cache/lj_rs_chain32.log
for v in "a2:0" "a2p1:1" "a2m1:-1" "a2p3:3.3" "a2m3:-3.3"; do
  tag=${v%%:*}; f=${v##*:}
  LJ_FEXT=$f python3 lj_rs_chain.py nl 32 3 4 2 $tag > cache/lj_rs_chain32_${tag}_nl.log 2>&1
  for c in keep drop; do BJ_C3H=$c python3 lj_rs_chain.py win 32 $tag > cache/lj_rs_chain32_${tag}_win_$c.log 2>&1; done
  grep -h "J=32: min" cache/lj_rs_chain32_${tag}_win_*.log | sed "s/^/$tag /" >> $L
done
echo done-d3 >> $L
