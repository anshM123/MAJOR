#!/bin/bash
# po_t15 (PARITY_ODD_CLOSE.md sect. 3): hc_t3.sh (copied) for ODD J = 15 with the new direct n3 nodes (row 400) in the z4 chain.
# usage: ./po_t15.sh variant...  -> cache/hc_chain_z4_po15_<variant>.*
cd "$(dirname "$0")"
J=15; ROWS="15:400"
EXM=1.02,1.05,1.0586,1.1,1.15,1.2,1.3,1.4,1.5,1.75,2,2.5,3,4,5
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ulimit -v 3000000
C="Q3_JEVEN=4 Q3_JODD=$J Q3_EXTRA_M=$EXM SM_CURV=row HC_NODES=0 Z4_JS=$J HC_Z3=$ROWS"
run() { env $C HC_TAG=po15_$1 "${@:2}" python3 hc_chain.py z4 300 > cache/hc_chain_z4_po15_$1.log 2>&1; }
for v in "$@"; do
  case $v in
    best) run best Z4_EST=best ;; dir) run dir Z4_EST=dir ;; cv) run cv Z4_EST=cv ;;
    zero) run zero Z4_EST=best Z4_EXTRAP=zero ;; lin) run lin Z4_EST=best Z4_EXTRAP=lin ;;
    convA) run convA Z4_DIRECT=0 Z4_NL=0 ;;
    jk*) run $v Z4_EST=best Z4_JK=${v#jk} ;;
  esac
done
# row injection (HC_INJ: both rows' direct nodes enter the odd K2 rows; z4's NL term does not act on odd J):
#   ./po_t15.sh inj inj_nodesp injjk0 ... injjk7
for v in "$@"; do
  case $v in
    inj) run inj HC_INJ=15:336,400 HC_Z3= ;; inj_nodesp) run inj_nodesp HC_INJ=15:336,400 HC_Z3= HC_DESPIKE_INJ=0 ;;
    injjk*) run $v HC_INJ=15:336,400 HC_Z3= HC_JK=${v#injjk} ;;
  esac
done
