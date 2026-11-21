#!/bin/bash
# H_OO (HYPOTHESES_CLOSE2.md): OO block at Delta_g = 300 in the z4 small-nu treatment -- variants + jackknife (oo_chain.py).
# usage: ./oo_t.sh PREFIX "J,J,..." [extra env]   -> cache/oo_chain_<PREFIX>_<variant>.{txt,pkl}
cd "$(dirname "$0")"
P=$1; JS=$2; shift 2
EXM=1.02,1.0395,1.0459,1.0523,1.0586,1.0649,1.0712,1.1,1.15,1.2,1.3,1.4,1.5,1.75,2,2.5,3,4,5
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
C="OO_J=$JS Q3_EXTRA_M=$EXM SM_CURV=row Z4_EST=best $*"
run() { env $C OO_TAG=${P}_$1 "${@:2}" python3 oo_chain.py > cache/oo_chain_${P}_$1.log 2>&1; }
run best OO_EST=best
run dir OO_EST=dir
run cv OO_EST=cv
run nodesp OO_EST=best OO_DESPIKE=0
run zero OO_EST=best Z4_EXTRAP=zero OO_RUN0=0
run lin OO_EST=best Z4_EXTRAP=lin OO_RUN0=0
run convA OO_DIRECT=0 Z4_DIRECT=0 Z4_NL=0
for c in -1000 -2000 -4000 -8000; do run c$c OO_EST=best Z4_EXTRAP=const Z4_CONST=$c OO_RUN0=0; done
for i in 0 1 2 3 4 5 6 7; do run jk$i OO_EST=best OO_JK=$i Z4_JK=$i; done
