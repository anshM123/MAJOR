#!/bin/bash
# po_t6 (PARITY_ODD_CLOSE.md sect. 3): oo_t.sh (copied) for J = 6 with the new direct 3x3 nodes (row 400): EE block (run0) and OO block.
# usage: ./po_t6.sh   -> cache/oo_chain_po6_<variant>.{txt,pkl,log}
cd "$(dirname "$0")"
EXM=1.0199,1.02,1.05,1.1,1.15,1.2,1.3,1.4,1.5,1.75,2,2.5,3,4,5
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ulimit -v 3000000
C="OO_J=6 Q3_EXTRA_M=$EXM SM_CURV=row Z4_EST=best OO_EESRC=oo OO_Z3EE=6:400"
run() { env $C OO_TAG=po6_$1 "${@:2}" python3 oo_chain.py > cache/oo_chain_po6_$1.log 2>&1; }
run best OO_EST=best
run dir OO_EST=dir
run cv OO_EST=cv
run zero OO_EST=best Z4_EXTRAP=zero OO_RUN0=0
run lin OO_EST=best Z4_EXTRAP=lin OO_RUN0=0
run convA OO_DIRECT=0 Z4_DIRECT=0 Z4_NL=0
for i in 0 1 2 3 4 5 6 7; do run jk$i OO_EST=best OO_JK=$i Z4_JK=$i; done
