#!/bin/bash
# q6: final analysis runs for S1_S4_FIXES.md (one process at a time).  Outputs in cache/q3_analysis_*.txt, cache/q5_tables_*.txt
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
EM=1.4,1.5,1.6,1.75,2.0,2.5
# main: K2, fixed C^_3, aligned rows, c_off model, all J
Q3_TAG=_main python3 q3_analysis.py 300 1000 > /dev/null 2>&1
Q3_TAG=_main Q3_EXTRA_M=$EM python3 q5_tables.py 300 10 > /dev/null 2>&1
Q3_TAG=_main Q3_EXTRA_M=$EM python3 q5_tables.py 1000 10 > /dev/null 2>&1
# l9 conventions for the off-diagonal (hold Gamma_off, no alignment) -- for comparison
Q3_TAG=_l9conv Q3_OFFSCALE=0 Q3_ALIGN=0 python3 q3_analysis.py 300 1000 > /dev/null 2>&1
# old C^_3 (scalar-g shift on k = 3)
Q3_TAG=_oldC3 Q3_C3SHIFT=scalar python3 q3_analysis.py 300 1000 > /dev/null 2>&1
# direct Delta = 400 rows (J = 16, 18, 20), RELMAX 0.05 (default) and 0.01 (x = 0.05 direct points dropped)
python3 q4_direct.py analyse 16 18 20 > /dev/null 2>&1
Q3_TAG=_direct Q3_DIRECT=1 python3 q3_analysis.py 300 > /dev/null 2>&1
Q3_TAG=_direct Q3_DIRECT=1 Q3_EXTRA_M=$EM python3 q5_tables.py 300 14 > /dev/null 2>&1
cp cache/q4_direct_grids.pkl cache/q4_direct_grids_relmax05.pkl; cp cache/q4_direct_compare.txt cache/q4_direct_compare_relmax05.txt
N1_RELMAX=0.01 python3 q4_direct.py analyse 16 18 20 > /dev/null 2>&1
cp cache/q4_direct_compare.txt cache/q4_direct_compare_relmax01.txt
Q3_TAG=_direct01 Q3_DIRECT=1 python3 q3_analysis.py 300 > /dev/null 2>&1
cp cache/q4_direct_grids_relmax05.pkl cache/q4_direct_grids.pkl; cp cache/q4_direct_compare_relmax05.txt cache/q4_direct_compare.txt
# small-x diagnostics: graviton correction switched off below x = 0.1 / 0.15
Q3_TAG=_xcut0.1 Q3_XCUT=0.1 python3 q3_analysis.py 300 1000 > /dev/null 2>&1
Q3_TAG=_xcut0.15 Q3_XCUT=0.15 python3 q3_analysis.py 300 1000 > /dev/null 2>&1
echo done
# J = 24, 28 (Delta = 1600 row only) with hybrid states at every twist-allowed m/Lam (crossing estimate, §2.6)
EX=1.02,1.03,1.04,1.05,1.06,1.07,1.08,1.09,1.1,1.12,1.15,1.2,1.3,1.4,1.5,1.6,1.75,2.0,2.5
Q3_TAG=_hiJ Q3_JEVEN=16,18,20,24,28 Q3_JODD=11,13 Q3_EXTRA_M=$EX python3 q5_tables.py 1000 16 > /dev/null 2>&1
Q3_TAG=_hiJ Q3_JEVEN=16,18,20,24,28 Q3_JODD=11,13 Q3_EXTRA_M=$EX python3 q5_tables.py 300 16 > /dev/null 2>&1
echo done hiJ
