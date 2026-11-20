#!/bin/bash
# after hc_queueA.sh: (o) finer edge nodes on (1040, 40), (1036, 36), (1031, 31) (the edge oscillation period is ~0.03 in x);
# (i) K2 beyond its criterion: direct bulk nodes on (336, 28), (336, 32) (n1 adaptive);
# (ii) J = 18 row 400 small-nu (n3); (iii) J = 4 / 5 small-nu spot nodes on row 400.  Resumable.
cd "$(dirname "$0")"
while pgrep -f "hc_queueA[.]sh" > /dev/null; do sleep 60; done
export OMP_NUM_THREADS=1
HC_METHOD=n1 HC_ADAPT=1 HC_M=13 HC_MREF=11 python3 hc_points.py 40 1040 0.8125,0.8375,0.8625
HC_METHOD=n1 HC_ADAPT=1 HC_M=13 HC_MREF=11 python3 hc_points.py 36 1036 0.8125,0.8375,0.8625
HC_METHOD=n1 HC_ADAPT=1 HC_M=12 HC_MREF=11 python3 hc_points.py 31 1031 0.8125,0.8375,0.8625
HC_METHOD=n1 HC_ADAPT=1 HC_M=12 HC_MREF=11 python3 hc_points.py 28 336 0.2,0.3,0.4,0.5,0.6
HC_METHOD=n1 HC_ADAPT=1 HC_M=12 HC_MREF=11 python3 hc_points.py 32 336 0.2,0.3,0.4,0.5,0.6
HC_METHOD=n3 HC_M=10 Z2_NUSPLIT=5 python3 hc_points.py 18 400 0.021,0.025,0.0293,0.0335,0.042,0.05,0.075,0.1
HC_METHOD=n3 HC_M=10 Z2_NUSPLIT=5 python3 hc_points.py 4 400 0.025,0.05,0.1
HC_METHOD=n3 HC_M=10 Z2_NUSPLIT=5 python3 hc_points.py 5 400 0.025,0.05,0.1
echo "QUEUE A2 DONE"
