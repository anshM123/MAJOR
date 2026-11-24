#!/bin/bash
# gaplarge (GAP_LARGESPIN.md): large-spin / large-impact-parameter pieces of Omega^EE positivity, uniform in Delta_g >= 1000.
# At most 2 processes; ~5 CPU-h in total.  Outputs in cache/ (git-ignored except .gitignore).
set -e
cd "$(dirname "$0")"
ulimit -v 3000000; export OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
python3 gl1_lrinf_uniform.py 2000 5000 > cache/gl1.log   # [certified] Lam^4-scaling form (envelope R^), direct checks at 2000/5000 (~10 s)
python3 gl2_lrinf_allgaps.py > cache/gl2.log             # [certified] LR-inf for EVERY Dg >= 1000 at beta0 = 0.37 / 0.245; uniform floor (~1 min)
python3 gl3_jstar.py > cache/gl3.log                     # [certified] J*(m; Dg) tables (Arb integer decisions)
for u in 4/9 100/169 1; do GL7_UMAX=$u python3 gl7_routea_uniform.py > cache/gl7_u$(echo $u|tr / _).log; done   # [D + certified] U3 (~3 min)
H13A_SIGNED=1 H13A_DGS=2000,5000,10000 H13A_S1=60,80,120 python3 gl5_routea_spot.py > cache/gl5_routea.log   # per-gap ROUTEA (~30 min)
GL4_DGS=1000,2000 python3 gl4_tol.py > cache/gl4.log; GL4_DGS=5000 python3 gl4_tol.py > cache/gl4_5000.log   # [DV] Theorem L' tolerances (~1 h)
python3 gl6_p4w.py > cache/gl6.log                       # [DV] Lemma P4W (uniform rate) on L' states at 1000/2000/5000
python3 gl8_thmW_chartB.py 1000 26,64,150,344 > cache/gl8_val.log   # check: reproduces agproof/cache/chi.json rows (~4 min)
python3 gl8_thmW_chartB.py 2000 > cache/gl8_2000.log &              # [DV] Theorem W at chart-B images (~20 min each)
python3 gl8_thmW_chartB.py 5000 > cache/gl8_5000.log; wait
python3 gl8_thmW_chartB.py 2000 edge > cache/gl8_2000_edge.log &    # twist-edge images (~10 min each)
python3 gl8_thmW_chartB.py 5000 edge > cache/gl8_5000_edge.log; wait
python3 gl4b_c2conv.py > cache/gl4_c2conv.log                 # [DV] node-convergence check of gl4 (~30 min)
