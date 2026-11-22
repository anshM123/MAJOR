#!/bin/bash
# cl_t1000 (RESULT_D_CLEANUP.md, H_edge at Delta_g = 1000, low J): the EE_NEARGAP main chain (hc_chain w2 1000 via cl_chain.py) with the
# direct edge nodes x = 0.8-0.9 of the rows ROWS injected (HC_IXMIN = 0.7; exempt from the Hampel rule), against the same chain without
# injection; twist-edge mass added.   usage: ./cl_t1000.sh J ROWS      -> cache/cl_chain_w2_k<J>_<variant>.*
cd "$(dirname "$0")"
J=$1; RE=$2
TW=$(python3 cl_tw.py 1000 $J)
EXM=$TW,1.01,1.02,1.05,1.1,1.15,1.2,1.3,1.4,1.5,1.75,2,2.5,3,4,5
if [ $((J % 2)) -eq 0 ]; then JJ="Q3_JEVEN=$J Q3_JODD=5"; else JJ="Q3_JEVEN=4 Q3_JODD=$J"; fi
C="$JJ Q3_EXTRA_M=$EXM SM_CURV=row"
run() { ./cl_pause.sh env $C HC_TAG=k${J}_$1 "${@:2}" python3 cl_chain.py hc w2 1000 > cache/cl_chain_w2_k${J}_$1.log 2>&1; }
run base
run edge HC_INJ=$J:$RE HC_IXMIN=0.7 HC_DESPIKE_INJ=0
run edge_desp HC_INJ=$J:$RE HC_IXMIN=0.7
S=$(python3 cl_shift.py w2 1000 $J k${J}_edge k${J}_base 2>/dev/null || echo 1); if python3 -c "import sys; sys.exit(0 if float('$S') > 0.003 else 1)"; then for i in 0 1 2 3 4 5 6 7; do run edgejk$i HC_INJ=$J:$RE HC_IXMIN=0.7 HC_DESPIKE_INJ=0 HC_JK=$i; done; fi   # edge jackknife only if the edge shifts a near-gap state by > 0.003 N_s
