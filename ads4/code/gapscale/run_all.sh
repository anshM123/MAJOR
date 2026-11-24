#!/bin/bash
# gapscale (GAP_SCALING.md): reproduce.  ~1.6 CPU-h in total, at most 2 processes; ulimit -v 3 GB, OMP_NUM_THREADS=1.
# Read-only imports: deshift/d_abel, spinfunc/ee1000_common (+ w2_chain, q3_analysis), lrtight/lrt_common + lrt3_rcert,
# neargap2/ng_chain + ng_bind (in-process; nothing is written there), spinfunc row engines (w1_rows, p3_odd_action, ...),
# agproof/ag_exact + p4win/p4w_common.  Every write goes to gapscale/cache.
set -e
cd "$(dirname "$0")"
ulimit -v 3000000; export OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
printf '*\n!.gitignore\n' > cache/.gitignore
python3 gs_weights.py 300 1000 2000 5000 > cache/weights.log            # §1   ~20 s
python3 gs_state.py > cache/state.log                                   # §2.1 ~20 s
python3 gs_lcscan.py > /dev/null                                        # §2.2 ~1 s   -> cache/lcscan.txt
for D in 1000 2000 5000; do python3 gs_lrinf.py $D > cache/lr_$D.log; done   # §3.1 ~10 s
python3 gs_lrball.py 1e6 > cache/lrball_part1.log 2>&1 || true          # §3.1 ~20 min; balls widen to x1.05 (stops with nan->assert only if wider)
#   (the published run was split: 998.5..1652 in lrball_part1.log, then `python3 gs_lrball.py 1e6 1550.04 > cache/lrball_resume.log`)
python3 gs_rows.py EE,C3H 5 2005 dense > cache/rows_J5_D2005.log         # §4   ~17 min (own-mass K2 row, Delta = 2005)
python3 gs_bind.py > cache/bind.txt                                     # §4   ~3 min
GS_PLO=0.002 GS_TAG=_full python3 gs_edge.py "1000:26,2000:52,500:13,1500:39,2000:26,5000:26,1000:48,2000:48" > cache/edge_full.log   # §3.3 ~6 min
GS_PLO=0.002 GS_TAG=_fixu python3 gs_edge.py "2000:26:1.0257,5000:26:1.0257,2000:48:1.04695" > /dev/null                          # §3.3 ~3 min
python3 gs_cover.py > /dev/null                                         # §5   <1 s   -> cache/cover.txt
