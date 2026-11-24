#!/bin/bash
# AG_PROOF.md: reproduce everything (~1.2 CPU-h; at most 2 processes).  Read-only imports from p4win/, p2con/, exactskel/, spinfunc/.
set -e
cd "$(dirname "$0")"
ulimit -v 3000000; export OMP_NUM_THREADS=1
printf '*\n!.gitignore\n' > cache/.gitignore
python3 ag_flat.py > cache/flat.out                               # §1-2  Lemma F checks; flat profile vs AG' (67 states)
python3 ag_unif.py > cache/unif.out                               # §1.1  uniform-rate / amplitude bounds of the flat waves
python3 ag_exact.py 40 48 64 100 > cache/exact.out                 # §3.1-3.2 exact scalar on vertical lines (complex nu)
AG_TRUNC=1 python3 ag_exact.py 40 100 >> cache/exact.out           #        truncation check (doubled M)
python3 ag_growth.py > cache/growth.out                              # §3.2  growth rates on vertical lines
python3 ag_pwmeasure.py > cache/pwmeasure.out                      # §7.1  beyond-edge plane-wave mass (MC)
AG_OFF=1.5 python3 ag_lock.py > cache/lock15.out & AG_OFF=3 python3 ag_lock.py > cache/lock3.out; wait   # §3.3 flat-locked fit (negative)
python3 ag_edge.py > cache/edge.out & python3 ag_edge2.py > cache/edge2.out; wait                       # §4, §5.3 edge contribution
python3 ag_chi.py > cache/chi.out & python3 ag_p4w.py > cache/p4w.out; wait                             # §5.2, §6
