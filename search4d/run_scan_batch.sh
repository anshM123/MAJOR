#!/bin/bash
# usage: run_scan_batch.sh dlt seed "kappas"  -> all t0 for this dlt, sequentially (each run < 20 min)
dlt=$1; seed=$2; kaps=${3:-"0.0 0.6 -0.6"}
for t0 in 0.0 0.3 0.6 1.0 1.5 2.5; do
 for kap in $kaps; do
  out=runs/v_d${dlt}_t${t0}_k${kap}_s${seed}.json
  [ -f $out ] && continue
  timeout 1150 python3 -I -u vaidya_scan.py $dlt $t0 $kap 8 $seed $out 2>/dev/null | grep "^{" >> runs/log_d${dlt}_s${seed}.txt
 done
done
echo DONE >> runs/log_d${dlt}_s${seed}.txt
