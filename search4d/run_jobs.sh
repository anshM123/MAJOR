#!/bin/bash
# usage: run_jobs.sh tag "dlt:t0:kappa:seed ..."   (each run < 20 min)
tag=$1; shift
for job in $@; do
  IFS=: read dlt t0 kap seed <<< "$job"
  out=runs/v_d${dlt}_t${t0}_k${kap}_s${seed}.json
  [ -f $out ] && continue
  timeout 1150 python3 -I -u vaidya_scan.py $dlt $t0 $kap 8 $seed $out 2>/dev/null | grep "^{" >> runs/log_$tag.txt
done
echo DONE >> runs/log_$tag.txt
