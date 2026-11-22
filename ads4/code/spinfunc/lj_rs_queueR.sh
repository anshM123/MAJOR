#!/bin/bash
# HLJ_CLOSE: reference (2e4, 32) depth-2 pieces (needed by lj_rs_norm for J = 32), then J = 32 normalisation, validation and window chain.
cd "$(dirname "$0")"
ulimit -v 3000000
export BJ_ADAPT=1 BJ_PRUNE=1e-12 BJ_KTOL=1e-6 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
while pgrep -f lj_rs_queue28.sh > /dev/null; do sleep 20; done
python3 lj_rs_eeq.py 20000 32 2 4 > cache/lj_rs_eeq_D20000_J32_d2.log 2>&1
python3 lj_rs_norm.py 32 336 2 5 1 > cache/lj_rs_norm_J32.log 2>&1
python3 lj_rs_val.py 32 2 5 1 > cache/lj_rs_val_J32.log 2>&1
source ./lj_rs_env32.sh
L=cache/lj_rs_chain32.log
for v in "a1:0" "a1p1:1" "a1m1:-1" "a1p3:3.3" "a1m3:-3.3"; do
  tag=${v%%:*}; f=${v##*:}
  LJ_FEXT=$f python3 lj_rs_chain.py nl 32 2 5 1 $tag > cache/lj_rs_chain32_${tag}_nl.log 2>&1
  for c in keep drop; do BJ_C3H=$c python3 lj_rs_chain.py win 32 $tag > cache/lj_rs_chain32_${tag}_win_$c.log 2>&1; done
  grep -h "J=32: min" cache/lj_rs_chain32_${tag}_win_*.log | sed "s/^/$tag /" >> $L
done
for g in -2000 -4000; do LJ_GSHIFT=$g LJ_TAGSFX=_g$g BJ_C3H=keep python3 lj_rs_chain.py win 32 a1 > cache/lj_rs_chain32_a1_g$g.log 2>&1
  grep -h "J=32: min" cache/lj_rs_chain32_a1_g$g.log | sed "s/^/a1 g$g /" >> $L; done
echo done >> $L
