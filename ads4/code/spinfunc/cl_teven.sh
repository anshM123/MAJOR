#!/bin/bash
# cl_teven (RESULT_D_CLEANUP.md): even J at Delta_g = 300 with the direct small-nu 3x3 nodes (row ROW, cl_oopts/oo_pts) --
#  (1) OO block + EE block (oo_chain via cl_chain.py oo; the oo_t.sh variant set)     -> cache/cl_oochain_cl<J>_<variant>.*
#  (2) EE block in hc_chain z4 with the same small-nu nodes (CL_Z3OO) and, optionally, the direct edge nodes x = 0.8-0.9 (HC_INJ,
#      HC_IXMIN = 0.7, exempt from the Hampel rule as in HYPOTHESES_CLOSE)           -> cache/cl_chain_z4_cl<J>_<variant>.*
# usage: ./cl_teven.sh J ROW [oo] [edge]       (runs the analysis while the edge-node worker is paused: cl_pause.sh)
cd "$(dirname "$0")"
J=$1; R=$2; shift 2
TW=$(python3 cl_tw.py 300 $J)
EXM=$TW,$(seq -s, 1.015 0.005 1.12),$(seq -s, 1.13 0.01 1.2),1.25,1.3,1.4,1.5,1.75,2,2.5,3,4,5     # dense near the gap (J >= 11 has no exact-skeleton states)
for part in "$@"; do
if [ "$part" == "oo" ]; then
  C="OO_J=$J Q3_EXTRA_M=$EXM SM_CURV=row Z4_EST=best OO_EESRC=oo OO_Z3EE=$J:$R"
  run() { ./cl_pause.sh env $C OO_TAG=cl${J}_$1 "${@:2}" python3 cl_chain.py oo > cache/cl_oochain_cl${J}_$1.log 2>&1; }
  run best OO_EST=best
  run dir OO_EST=dir
  run cv OO_EST=cv
  run nodesp OO_EST=best OO_DESPIKE=0
  run zero OO_EST=best Z4_EXTRAP=zero OO_RUN0=0
  run lin OO_EST=best Z4_EXTRAP=lin OO_RUN0=0
  run convA OO_DIRECT=0 Z4_DIRECT=0 Z4_NL=0
  for c in -1000 -2000 -4000 -8000; do run c$c OO_EST=best Z4_EXTRAP=const Z4_CONST=$c OO_RUN0=0; done
  for i in 0 1 2 3 4 5 6 7; do run jk$i OO_EST=best OO_JK=$i Z4_JK=$i; done
fi
if [ "$part" == "edge" ]; then
  Z3OO=$J:$R; [ "$CL_Z3SRC" == "hc" ] && Z3OO=      # CL_Z3SRC=hc: small-nu EE nodes from hc_pts (J = 18-22), else the 3x3 oo points
  C="Q3_JEVEN=$J Q3_JODD=5 Q3_EXTRA_M=$EXM SM_CURV=row Z4_EST=best Z4_JS=$J HC_Z3=$J:$R CL_Z3OO=$Z3OO"
  run() { ./cl_pause.sh env $C HC_TAG=cl${J}_$1 "${@:2}" python3 cl_chain.py hc z4 300 > cache/cl_chain_z4_cl${J}_$1.log 2>&1; }
  run ebase
  run edge HC_INJ=$J:$R HC_IXMIN=0.7 HC_DESPIKE_INJ=0
  run edge_desp HC_INJ=$J:$R HC_IXMIN=0.7
  S=$(python3 cl_shift.py z4 300 $J cl${J}_edge cl${J}_ebase 2>/dev/null || echo 1); if python3 -c "import sys; sys.exit(0 if float('$S') > 0.003 else 1)"; then for i in 0 1 2 3 4 5 6 7; do run edgejk$i HC_INJ=$J:$R HC_IXMIN=0.7 HC_DESPIKE_INJ=0 HC_JK=$i; done; fi   # edge jackknife only if the edge shifts a near-gap state by > 0.003 N_s
fi
done
