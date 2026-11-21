#!/bin/bash
# Task 3 (Delta_g = 300 small-nu, z4 protocol) analysis for one even J with hc small-nu rows: variants + jackknife.
# usage: ./hc_t3.sh J "J:D[,D]" [extra env]   -> cache/hc_chain_z4_t3_<J>_<variant>.*
cd "$(dirname "$0")"
J=$1; ROWS=$2; shift 2
EXM=1.02,1.05,1.058,1.0656,1.07,1.072,1.1,1.15,1.2,1.3,1.4,1.5,1.75,2,2.5,3,4,5
C="Q3_JEVEN=$J Q3_JODD=5 Q3_EXTRA_M=$EXM SM_CURV=row HC_NODES=0 Z4_JS=$J HC_Z3=$ROWS $*"
run() { env $C HC_TAG=t3_${J}_$1 "${@:2}" python3 hc_chain.py z4 300 > cache/hc_chain_z4_t3_${J}_$1.log 2>&1; }
run best Z4_EST=best
run dir Z4_EST=dir
run cv Z4_EST=cv
run zero Z4_EST=best Z4_EXTRAP=zero
run lin Z4_EST=best Z4_EXTRAP=lin
run convA Z4_DIRECT=0 Z4_NL=0
for c in -1000 -2000 -4000 -8000; do run c$c Z4_EST=best Z4_EXTRAP=const Z4_CONST=$c; done
for i in 0 1 2 3 4 5 6 7; do run jk$i Z4_EST=best Z4_JK=$i; done
