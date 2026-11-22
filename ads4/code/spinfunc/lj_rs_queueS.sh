#!/bin/bash
# HLJ_CLOSE §4: structure-part stress of the (336, 28/32) windows (the unexplained structure offsets of §3, scaled), after lane 2.
cd "$(dirname "$0")"
ulimit -v 3000000
while [ ! -f cache/lj_rs_queueG.done ]; do sleep 30; done
L=cache/lj_rs_chainS.log
( source ./lj_rs_env28.sh
  for v in "s1:-50,-60,-25" "s3:-150,-180,-75" "s3p:150,180,75" "s10:-500,-600,-250"; do
    tag=${v%%:*}; d=${v##*:}
    LJ_SSTRUCT=$d python3 lj_rs_chain.py nl 28 3 4 2 $tag > cache/lj_rs_chain28_${tag}_nl.log 2>&1
    BJ_C3H=keep python3 lj_rs_chain.py win 28 $tag > cache/lj_rs_chain28_${tag}_win_keep.log 2>&1
    grep -h "J=28: min" cache/lj_rs_chain28_${tag}_win_keep.log | sed "s/^/J28 $tag ($d) /" >> $L
  done )
( source ./lj_rs_env32.sh
  for v in "s1:-224,-83,-82" "s3:-672,-249,-246" "s3p:672,249,246"; do
    tag=${v%%:*}; d=${v##*:}
    LJ_SSTRUCT=$d python3 lj_rs_chain.py nl 32 2 5 1 $tag > cache/lj_rs_chain32_${tag}_nl.log 2>&1
    BJ_C3H=keep python3 lj_rs_chain.py win 32 $tag > cache/lj_rs_chain32_${tag}_win_keep.log 2>&1
    grep -h "J=32: min" cache/lj_rs_chain32_${tag}_win_keep.log | sed "s/^/J32 $tag ($d) /" >> $L
  done )
echo done >> $L
# window thimble nodes (n4, nu = 5.0, 6.67) included in the J = 28 validation
LJ_VMIN=4.0 LJ_VTAG=_win OMP_NUM_THREADS=1 python3 lj_rs_val.py 28 3 4 2 > cache/lj_rs_val_J28_win.log 2>&1
echo done-val >> cache/lj_rs_chainS.log
