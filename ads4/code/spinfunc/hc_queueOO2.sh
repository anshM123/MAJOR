#!/bin/bash
# H_OO at larger J (comparative: OO vs EE ++ in the same K2 protocol): OO rows J = 18 at Delta = 336, 400 (11 nodes);
# waits for hc_queueB.sh.  Resumable.
cd "$(dirname "$0")"
while pgrep -f "hc_queueB[.]sh" > /dev/null || pgrep -f "hc_queueOO[.]sh" > /dev/null; do sleep 60; done
export OMP_NUM_THREADS=1 HC_XS=0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8 HC_DS=336,400
python3 hc_oo_rows.py EE 18; python3 hc_oo_rows.py C3H 18
echo "QUEUE OO2 DONE"
