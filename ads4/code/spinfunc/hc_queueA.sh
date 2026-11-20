#!/bin/bash
# H_edge: direct (n1 adaptive, shared seeds, Q_J rescaling, Fourier floor) nodes on the light-cone-edge window x in (0.8, 0.9] of the
# rows that carry the near-gap states (resumable: finished nodes are skipped; re-run to resume after a restart)
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 HC_METHOD=n1 HC_ADAPT=1 HC_MREF=11
HC_M=13 python3 hc_points.py 40 1040 0.8,0.825,0.85,0.875,0.9
HC_M=12 python3 hc_points.py 31 1031 0.8,0.825,0.85,0.875,0.9
HC_M=13 python3 hc_points.py 36 1036 0.8,0.825,0.85,0.875,0.9
HC_M=12 python3 hc_points.py 24 336 0.8,0.825,0.85,0.875,0.9
HC_M=12 python3 hc_points.py 17 1017 0.8,0.825,0.85,0.875,0.9
HC_M=12 python3 hc_points.py 16 400 0.8,0.825,0.85,0.875,0.9
for D in 200 400 800 1600; do HC_M=10 python3 hc_points.py 5 $D 0.8,0.85,0.9,0.95; done
echo "QUEUE A DONE"
