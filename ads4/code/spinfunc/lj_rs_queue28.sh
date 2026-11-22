#!/bin/bash
# HLJ_CLOSE §4: (336, 28) window with the analytically normalised engine (a^2, depth 3, N = 4), production + stress variants.
cd "$(dirname "$0")"
ulimit -v 3000000
source ./lj_rs_env28.sh
L=cache/lj_rs_chain28.log
for v in "a2:0" "a2p1:1" "a2m1:-1" "a2p3:3.3" "a2m3:-3.3"; do
  tag=${v%%:*}; f=${v##*:}
  LJ_FEXT=$f python3 lj_rs_chain.py nl 28 3 4 2 $tag > cache/lj_rs_chain28_${tag}_nl.log 2>&1
  for c in keep drop; do BJ_C3H=$c python3 lj_rs_chain.py win 28 $tag > cache/lj_rs_chain28_${tag}_win_$c.log 2>&1; done
  grep -h "J=28: min" cache/lj_rs_chain28_${tag}_win_*.log | sed "s/^/$tag /" >> $L
done
# a^1 only (depth 3, Kr = 1) for the order-by-order table
python3 lj_rs_chain.py nl 28 3 4 1 a1 > cache/lj_rs_chain28_a1_nl.log 2>&1
for c in keep drop; do BJ_C3H=$c python3 lj_rs_chain.py win 28 a1 > cache/lj_rs_chain28_a1_win_$c.log 2>&1; done
grep -h "J=28: min" cache/lj_rs_chain28_a1_win_*.log | sed "s/^/a1 /" >> $L
# failure scale: uniform Gamma_c shifts on top of the production engine
for g in -2000 -4000 -8000; do LJ_GSHIFT=$g LJ_TAGSFX=_g$g BJ_C3H=keep python3 lj_rs_chain.py win 28 a2 > cache/lj_rs_chain28_a2_g$g.log 2>&1
  grep -h "J=28: min" cache/lj_rs_chain28_a2_g$g.log | sed "s/^/a2 g$g /" >> $L; done
echo done >> $L
