#!/bin/bash
# LARGESPIN_ANALYTIC.md: reproduce every table (<= 2 processes; ~1.5 CPU-h total).
set -e
cd "$(dirname "$0")"
ulimit -v 3000000; export OMP_NUM_THREADS=1
python3 t_core.py  > cache/t_core.txt         # OmLR (Mehler-Dirichlet) vs mpmath hyp2f1
python3 t_hyb.py   > cache/t_hyb.txt          # hyb_fast vs spinfunc/ee1000_common.hyb
python3 ls_barrier.py 1000                    # OmLR_2(t), OmLR_3(t) table, monotonicity        -> cache/barrier_1000.txt
python3 ls_struct.py 1000 > cache/struct_run.log   # tolerance profile + structures A-G      -> cache/struct_1000.{txt,json}
python3 ls_xterm.py 1000  > cache/xterm_run.log    # de-shift (X - psi2) term                  -> cache/xterm_1000.txt
python3 ls_curvfloor.py 1000 > cache/curvfloor_run.log   # floor with computed curvature       -> cache/curvfloor_1000.txt
python3 ls_curvbound.py 1000 > cache/curvbound_run.log   # Lemma C' (hybrid floor sign-free)  -> cache/curvbound_1000.txt
python3 ls_strip.py 1000                      # H_bulk relative form, strips                       -> cache/strip_1000.txt
python3 t_xsigned.py > cache/t_xsigned.txt    # X does not act on the stripped ring (sanity)
python3 ls_agscan.py 1000 2:50,2:150,3:150,5:250,1.5:225 > cache/agscan_run.log   # AG on H_R' states -> cache/agscan_1000.txt
