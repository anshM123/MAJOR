#!/bin/bash
# gaplarge (GAP_LARGESPIN.md): large-spin / large-impact-parameter pieces of Omega^EE positivity, uniform in Delta_g >= 1000.
# At most 2 processes; ~1 CPU-h in total.  Outputs in cache/ (git-ignored except .gitignore).
set -e
cd "$(dirname "$0")"
ulimit -v 3000000; export OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
python3 gl1_lrinf_uniform.py 2000 5000 > cache/gl1.log   # [certified] Lam^4-scaling form (envelope R^), direct checks at 2000/5000 (~10 s)
python3 gl2_lrinf_allgaps.py > cache/gl2.log             # [certified] LR-inf for EVERY Dg >= 1000 at beta0 = 0.37 / 0.245; uniform floor (~1 min)
python3 gl3_jstar.py > cache/gl3.log                     # [certified] J*(m; Dg) tables (Arb integer decisions)
python3 gl7_routea_uniform.py > cache/gl7.log            # [D + certified] ROUTEA window uniform in Dg >= 1000 (~1 min)
H13A_SIGNED=1 H13A_DGS=2000,5000,10000 H13A_S1=60,80,120 python3 gl5_routea_spot.py > cache/gl5_routea.log   # per-gap ROUTEA (~30 min)
python3 gl4_tol.py > cache/gl4.log                       # [DV] Theorem L' tolerances at 1000/2000/5000 (~40 min)
python3 gl6_p4w.py > cache/gl6.log                       # [DV] Lemma P4W (uniform rate) on L' states at 1000/2000/5000
