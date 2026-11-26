#!/bin/bash
# tailaudit/ta_rerun.sh: re-run the TAIL_CERT scripts in a private mirror (cache/mirror/ads4), so nothing in code/tailcert/cache is
# overwritten. Read-only inputs (contour, epswcert, lrtight, prl) are symlinked. Outputs: cache/mirror/ads4/code/tailcert/cache/.
set -e
ulimit -v 4000000; export OMP_NUM_THREADS=1
HERE=$(cd "$(dirname "$0")" && pwd); ADS4=$(cd "$HERE/../.." && pwd)
M=$HERE/cache/mirror/ads4
mkdir -p $M/code/tailcert/cache
for d in contour epswcert lrtight; do ln -sfn $ADS4/code/$d $M/code/$d; done
ln -sfn $ADS4/prl $M/prl
cp $ADS4/code/tailcert/*.py $M/code/tailcert/
cd $M/code/tailcert
python3 -B tc_exact.py > cache/run_exact.log 2>&1
python3 -B tc_power.py > cache/run_power.log 2>&1
python3 -B tc_cutoff.py > cache/run_cutoff.log 2>&1
python3 -B tc_lemmaR.py > cache/run_lemmaR.log 2>&1
( python3 -B tc_assemble.py 20 200 > cache/run_assemble_20.log 2>&1; python3 -B tc_k3.py 20 200 > cache/run_k3_20.log 2>&1 ) &
( python3 -B tc_remainder.py 20 200 > cache/run_rem_20.log 2>&1; python3 -B tc_remainder.py 18 200 > cache/run_rem_18.log 2>&1 ) &
( export TC_L0TAB=18; python3 -B tc_assemble.py 18 200 > cache/run_assemble_18.log 2>&1; python3 -B tc_k3.py 18 200 > cache/run_k3_18.log 2>&1 ) &
wait
echo done
