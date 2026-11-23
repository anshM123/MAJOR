#!/bin/bash
# LARGESPIN_L.md: reproduce everything (< 0.5 CPU-h, one process at a time). Outputs in cache/.
ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$(dirname "$0")"
python3 lsl_frame.py        # §2.1  gauge identity, Laplace/Euler/Riesz integrals        -> cache/frame_checks.txt   (~20 min, mpmath)
python3 lsl_lrinf.py        # §3    Lemma LR-inf: identities, remainder bound, direct check -> cache/lrinf.txt
python3 lsl_lrinf_grid.py   # §3    fine beta grid of the hybrid-level remainder           -> cache/lrinf_grid.txt
python3 lsl_hybcmp.py       # §3.3  hybrid vs F_LR; J*(m)                                  -> cache/hybcmp.txt
python3 lsl_smear.py        # §2.6, §4.1 Bernstein smearing vs CHMRS; lightcone uniformity -> cache/smear.txt
python3 lsl_window.py       # §4.2  window tolerance map, de-shift vs raw, J = 26..1e6      -> cache/window_map.{txt,json}
python3 lsl_ringtol.py      # §4.3  ring / AdS-local tolerances, large-beta profile         -> cache/ringtol.txt
