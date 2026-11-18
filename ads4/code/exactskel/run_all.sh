#!/bin/bash
# Exact scalar-skeleton heavy moments on the near-gap states (EXACT_SKELETON_NEARGAP.md).  <= 2 cores.
# Needs (read-only): ../omegau (common.py, adsaction.py, remk.py, u4_scan.py, cache/rcoef_Df*.pkl), ../robust, ../uniform (skel.py, caches).
# Usage: bash run_all.sh [Dg]        (default Dg = 100; about 50 min per skeleton on one core)
set -e
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
DG=${1:-100}
mkdir -p cache
# validation: Mack polynomials; Mellin moments vs the CHMRS order-1/m^8 series at small nu/m; contour/parameter independence at large nu/m
python3 mack.py              > cache/v0_mack.txt
python3 t_conv.py 120.3 0 3  > cache/v1_conv.txt
python3 t_check2.py          > cache/v2_series.txt
python3 t_largenu.py         > cache/v3_largenu.txt
# exact moments on the states (two skeletons in parallel), then Df = 2.5
python3 e1_states.py $DG 3 > cache/e1_${DG}_Df3.log 2>&1 &
python3 e1_states.py $DG 2 > cache/e1_${DG}_Df2.log 2>&1
wait
python3 e1_states.py $DG 2.5 > cache/e1_${DG}_Df2.5.log 2>&1 &
python3 e2_errcheck.py $DG 3 > /dev/null
python3 e4_quad.py $DG 3 > /dev/null
wait
python3 e3_action.py $DG 3 2 2.5 > /dev/null
python3 e5_compare.py $DG > /dev/null
# --- second/third pass: high-J, Arb method for Delta_g = 300/1000, checks, NF02, bounds (see EXACT_SKELETON_NEARGAP.md section 8)
MS=1.01,1.02,1.03,1.04,1.05,1.06,1.07,1.08,1.09,1.1,1.12,1.15,1.2,1.3
E1_PKL=_hiJ E1_JLIST=11,12,13,14,16,18,20,24,28 E1_MSET=1.11,1.12,1.15,1.2,1.25,1.3,1.35,1.4,1.5,1.6,1.75,2.0,2.25,2.5,3.0 python3 e1_states.py 100 3 3.0 > cache/e1_100_Df3_hiJ.log 2>&1 &
python3 e8_states_arb.py 300 3 2.0 10 > cache/e1_300_Df3.log 2>&1; wait
python3 e8_states_arb.py 300 2 2.0 10 > cache/e1_300_Df2.log 2>&1 &
python3 e8_states_arb.py 1000 3 1.5 10 > cache/e1_1000_Df3.log 2>&1; wait
python3 e8_states_arb.py 1000 2 1.5 10 > cache/e1_1000_Df2.log 2>&1 &
for D in 1000 300; do E1_PKL=_hiJ E1_JLIST=11,12,13,14,16,18,20 E1_MSET=$MS python3 e8_states_arb.py $D 3 1.5 > cache/e1_${D}_Df3_hiJ.log 2>&1; done; wait
python3 e2b_arbcheck.py 100 3 1.5 > /dev/null & python3 e2b_arbcheck.py 100 2 1.5 > /dev/null; wait
for D in 100 300 1000; do python3 e3_action.py $D 3 2 > /dev/null; python3 e7_nf.py $D > /dev/null; done
python3 e6_bound.py 100 1.08 1.09; python3 e6_bound.py 300 1.05 1.06; python3 e6_bound.py 1000 1.01
for D in 100 300 1000; do python3 e9_bound_nf.py $D 1.0; done
# --- section 9: even kernel S_ff vs K* (normalisations, thresholds with finite-m <TTTT> corrections, bounds)
python3 evenkernel.py > /dev/null
for D in 100 300 1000; do E_KERNEL=Kstar python3 e10_kernel.py $D > /dev/null; done; E_KERNEL=Sff python3 e10_kernel.py 100 > /dev/null
E_KERNEL=Kstar python3 e6_bound.py 100 1.08; E_KERNEL=Kstar python3 e6_bound.py 300 1.05; E_KERNEL=Kstar python3 e6_bound.py 1000 1.01
for D in 100 300 1000; do E_KERNEL=Kstar python3 e9_bound_nf.py $D 1.0; done
