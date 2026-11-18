"""
p3_run: production driver for p3_odd_action (resumable; one process per (kind, J-list) so that cache files are not shared).
For each J, Delta row and x in XS: the state point (nu = x m(Delta, J)) and the heavy reference (Delta_ref = 2e4, same nu, same kind).
Usage: OMP_NUM_THREADS=1 P3_DS=100,200,400,800,1600,3200 P3_XS=... python3 p3_run.py kind J [J ...]
"""
import os, sys
sys.dont_write_bytecode = True
os.environ.setdefault('L1_DEG', '6')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import p3_odd_action as P

XS = [float(v) for v in os.environ.get('P3_XS', '0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8').split(',')]
DS = [float(v) for v in os.environ.get('P3_DS', '1600,3200,800,400,200,100').split(',')]
DREF = float(os.environ.get('P3_DREF', '20000'))

if __name__ == '__main__':
    kind = sys.argv[1]
    for J in [int(a) for a in sys.argv[2:]]:
        for D in DS:
            m = P.mphys(D, J)
            for x in XS:
                P.get(kind, J, D, [x*m])
                P.get(kind, J, DREF, [x*m])
        print(f"done {kind} J={J}", flush=True)
