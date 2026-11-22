#!/bin/bash
# cl_todd (RESULT_D_CLEANUP.md): odd J at Delta_g = 300 -- direct small-nu nodes (x <= 0.2) and direct edge nodes (x >= 0.7) by row
# injection into the odd K2 rows (HC_INJ; z4's non-local term does not act on odd J), hc_chain z4 via cl_chain.py.
# usage: ./cl_todd.sh J "SMALLROWS" "EDGEROWS" [small] [edge]     -> cache/cl_chain_z4_cl<J>_<variant>.*
cd "$(dirname "$0")"
J=$1; RS=$2; RE=$3; shift 3
TW=$(python3 cl_tw.py 300 $J)
EXM=$TW,1.02,1.05,1.1,1.15,1.2,1.3,1.4,1.5,1.75,2,2.5,3,4,5
C="Q3_JEVEN=4 Q3_JODD=$J Q3_EXTRA_M=$EXM SM_CURV=row Z4_JS=$J Z4_EST=best"
run() { ./cl_pause.sh env $C HC_TAG=cl${J}_$1 "${@:2}" python3 cl_chain.py hc z4 300 > cache/cl_chain_z4_cl${J}_$1.log 2>&1; }
for part in "$@"; do
if [ "$part" == "small" ]; then
  run convA
  run inj HC_INJ=$J:$RS HC_IXMAX=0.2
  run inj_nodesp HC_INJ=$J:$RS HC_IXMAX=0.2 HC_DESPIKE_INJ=0
  run inj_drop HC_INJ=$J:$RS HC_IXMIN=0.0225 HC_IXMAX=0.2
  for i in 0 1 2 3 4 5 6 7; do run injjk$i HC_INJ=$J:$RS HC_IXMAX=0.2 HC_JK=$i; done
fi
if [ "$part" == "edge" ]; then
  # small-nu (where present) + edge rows; HC_INJ takes all direct nodes of the listed rows; edge nodes exempt from the Hampel rule
  run edge HC_INJ=$J:$RE HC_DESPIKE_INJ=0
  run edge_desp HC_INJ=$J:$RE
  run edgeonly HC_INJ=$J:$RE HC_IXMIN=0.7 HC_DESPIKE_INJ=0
  for i in 0 1 2 3 4 5 6 7; do run edgejk$i HC_INJ=$J:$RE HC_DESPIKE_INJ=0 HC_JK=$i; done
fi
done
