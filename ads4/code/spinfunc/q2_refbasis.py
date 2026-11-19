"""q2_refbasis (S1_S4_FIXES.md §2): C^_2 points at the heavy reference itself, x_ref = nu/m_ref = 0.3 ... 0.8, needed by l8_analysis.basis
(real helicity basis at Delta_ref from the order-0 action at x >= 0.15) for J >= 14 (J <= 12 have them from l8).  Cache as q2_run (tag _c).
Usage: OMP_NUM_THREADS=1 python3 q2_refbasis.py 14,16,18,20"""
import os, sys
os.environ['L1_DEG'] = '6'; os.environ['Q2_TAG'] = '_c'
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import q2_run as Q
P = Q.setup_even()
for J in [int(a) for a in sys.argv[1].split(',')]:
    mr = Q.mphys(Q.DREF, J)
    P.get('C2', J, Q.DREF, [x*mr for x in (0.3, 0.4, 0.5, 0.6, 0.7, 0.8)], verbose=False)
    print(f"done ref basis J={J}", flush=True)
