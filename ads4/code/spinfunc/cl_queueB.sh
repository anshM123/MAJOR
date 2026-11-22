#!/bin/bash
# cl_queueB (RESULT_D_CLEANUP.md, task 2): H_edge -- direct n1 nodes (adaptive proposal, shared seeds, Q_J rescaling, Fourier floor;
# hc_points n1 path unchanged) on the light-cone-edge window x = 0.8 ... 0.9 of the rows that carry the near-gap states.
# Delta_g = 300 first (rows 400 / 336 / 200), then the low-J rows that carry the Delta_g = 1000 gap states (1600 / 800 / 400).
# Even J: w0 2x2 EE integrand + C3H; odd J: 1x1 odd integrand.  -> cache/cl_pts_<kind>_J<J>_D<D>.pkl.  Resumable.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CL_MODE=hc HC_METHOD=n1 HC_ADAPT=1 HC_MREF=11 HC_M=11
ulimit -v 3000000
XE=0.8,0.825,0.85,0.875,0.9
e() { python3 cl_points.py $1 $2 $XE > cache/cl_qB_J$1_$2.log 2>&1; }
# Delta_g = 300 (lowest-margin first)
e 4 400; e 6 400; e 8 400; e 9 400; e 9 200; e 10 400; e 11 400; e 12 400; e 15 400; e 15 336; e 14 400; e 7 400; e 7 200; e 13 400
e 18 336; e 18 400; e 20 400; e 22 336
# Delta_g = 1000 low J (rows bracketing m ~ 1000)
e 4 1600; e 6 1600; e 7 1600; e 7 800; e 8 1600; e 9 1600; e 9 800; e 10 1600; e 11 800; e 13 800; e 12 1600; e 12 800
e 14 1600; e 14 800; e 16 1600; e 16 800
echo "QUEUE CL-B DONE"
