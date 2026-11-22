#!/bin/bash
# cl_driver5 (RESULT_D_CLEANUP.md): odd J = 9, 5, 7 at Delta_g = 300 with the direct small-nu nodes on BOTH rows that bracket the near-gap
# states (200 and 400; cl_queueA2 adds the 200 row), plus the edge nodes.  Tags cl<J>r2_* (J = 9, 7: cl_todd.sh variants; J = 5: edge
# x <= 0.92 as HYPOTHESES_CLOSE §2.4).
cd "$(dirname "$0")"
wA() { until grep -q "^# done" cache/cl_qA_J$1_$2.log 2>/dev/null; do sleep 30; done; }
while pgrep -f "cl_driver[.]sh" > /dev/null; do sleep 30; done
EXM5=1.0166,$(seq -s, 1.015 0.005 1.12),$(seq -s, 1.13 0.01 1.2),1.25,1.3,1.4,1.5,1.75,2,2.5,3,4,5
for J in 9 5 7; do
  wA $J 200
  if [ $J -eq 5 ]; then
    C="Q3_JEVEN=4 Q3_JODD=5 Q3_EXTRA_M=$EXM5 SM_CURV=row Z4_JS=5 Z4_EST=best"
    run() { ./cl_pause.sh env $C HC_TAG=cl5r2_$1 "${@:2}" python3 cl_chain.py hc z4 300 > cache/cl_chain_z4_cl5r2_$1.log 2>&1; }
    run convA
    run inj HC_INJ=5:200,400 HC_IXMAX=0.2
    run inj_nodesp HC_INJ=5:200,400 HC_IXMAX=0.2 HC_DESPIKE_INJ=0
    for i in 0 1 2 3 4 5 6 7; do run injjk$i HC_INJ=5:200,400 HC_IXMAX=0.2 HC_JK=$i; done
    run edge HC_INJ=5:200,400,800,1600 HC_IXMAX=0.92 HC_DESPIKE_INJ=0
  else
    sed "s/HC_TAG=cl\${J}_/HC_TAG=cl\${J}r2_/; s/cl_chain_z4_cl\${J}_/cl_chain_z4_cl\${J}r2_/; s/cl\${J}_edge cl\${J}_inj_nodesp/cl\${J}r2_edge cl\${J}r2_inj_nodesp/" cl_todd.sh > cl_todd_r2.sh
    chmod +x cl_todd_r2.sh
    ./cl_todd_r2.sh $J 200,400 200,400 small edge
  fi
done
echo "CL-DRIVER5 DONE"
