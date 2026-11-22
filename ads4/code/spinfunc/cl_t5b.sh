#!/bin/bash
# cl_t5b (RESULT_D_CLEANUP.md): J = 5 at Delta_g = 300, edge nodes restricted to x <= 0.92 (x = 0.95 has ~100 % RQMC error on the
# 200/400 rows; HYPOTHESES_CLOSE §2.4), small-nu spot nodes of row 400 included.   -> cache/cl_chain_z4_cl5_e92*.*
cd "$(dirname "$0")"
J=5; EXM=1.0166,1.02,1.05,1.1,1.15,1.2,1.3,1.4,1.5,1.75,2,2.5,3,4,5
C="Q3_JEVEN=4 Q3_JODD=5 Q3_EXTRA_M=$EXM SM_CURV=row Z4_JS=5 Z4_EST=best HC_INJ=5:200,400,800,1600 HC_IXMAX=0.92"
run() { ./cl_pause.sh env $C HC_TAG=cl5_$1 "${@:2}" python3 cl_chain.py hc z4 300 > cache/cl_chain_z4_cl5_$1.log 2>&1; }
run e92 HC_DESPIKE_INJ=0
run e92_desp
for i in 0 1 2 3 4 5 6 7; do run e92jk$i HC_DESPIKE_INJ=0 HC_JK=$i; done
