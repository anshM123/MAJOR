#!/bin/bash
# z6: replica-jackknife propagation of the z2 node errors through the z4 chain (8 leave-one-replica-out runs per estimator)
# usage: ./z6_jk.sh EST EXTRAP TAGPREFIX   -> cache/z4_chain_<TAGPREFIX>_jk<i>.txt ; summary by z6 lines in stdout
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 Z4_JS=24,28 Q3_JEVEN=24,28 Q3_JODD=5 Q3_EXTRA_M=1.02,1.05,1.1,1.15,1.2,1.3,1.4,1.5,1.75,2,2.5,3,4,5
for i in 0 1 2 3 4 5 6 7; do
  Z4_EST=$1 Z4_EXTRAP=$2 Z4_JK=$i Z4_TAG=${3}_jk$i python3 z4_chain.py 300 2>&1 | grep "J=2[48]:" | sed "s/^/jk$i /"
done
