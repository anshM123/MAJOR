#!/bin/bash
# nu < 8.3 window at Delta_g = 300, J = 28, 32 (HYPOTHESES_CLOSE2.md §3): hc_chain z4 with the direct bulk sub-rows (as HYPOTHESES_CLOSE
# §5), thimble node(s) at nu = 5.0 from a chosen source (Z3_PREFER), and a uniform Gamma_c := c below the lowest node (Z4_EXTRAP=const).
# usage: ./oo_window.sh J SRC "c1 c2 ..." [NUMIN]   (SRC = z3 prefer list, e.g. x | x,s2 | x,oom; NUMIN = 8 drops all thimble nodes)
cd "$(dirname "$0")"
J=$1; SRC=$2; CS=$3; NUMIN=${4:-0}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
if [ "$J" = 28 ]; then ROWS=336,500,1600; EXM=1.091,1.1,1.15,1.2,1.3,1.4,1.5,1.75,2,2.5,3,4,5; else ROWS=336,1040; EXM=1.102,1.1,1.15,1.2,1.3,1.4,1.5,1.75,2,2.5,3,4,5; fi
B="HC_INJ=$J:336 HC_IXMIN=0.12 HC_IXMAX=0.81 HC_DESPIKE_INJ=0 HC_NODES=0 Z4_EST=best Z4_JS=$J Q3_JEVEN=$J Q3_JODD=5 W2_ROWS=$ROWS Q3_EXTRA_M=$EXM SM_CURV=row Z3_PREFER=$SRC Z4_NUMIN=$NUMIN"
S=$(echo $SRC | tr ',' '_')
env $B HC_TAG=oowin_${J}_${S}_n${NUMIN}_hold python3 hc_chain.py z4 300 > /dev/null 2>&1
for c in $CS; do env $B Z4_EXTRAP=const Z4_CONST=$c HC_TAG=oowin_${J}_${S}_n${NUMIN}_c$c python3 hc_chain.py z4 300 > /dev/null 2>&1; done
if [ -n "$JK" ]; then for i in 0 1 2 3 4 5 6 7; do env $B Z4_JK=$i HC_TAG=oowin_${J}_${S}_n${NUMIN}_jk$i python3 hc_chain.py z4 300 > /dev/null 2>&1; done; fi
echo "WINDOW $J $SRC done"
