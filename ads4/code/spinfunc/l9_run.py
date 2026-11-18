"""
l9_run (doc 23): driver for l9_points: for each Delta in L9_DS and x in L9_XS, the state point (Delta, nu = x m) and the heavy reference
(Delta_ref = L9_DREF, same nu), kind and J from argv.  Resumable (l9_points cache).
Usage: L1_DEG=8 L9_DS=100,150 python3 l9_run.py kind J [J ...]
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import l9_points as P

XS = [float(v) for v in os.environ.get('L9_XS', '0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8').split(',')]
DS = [float(v) for v in os.environ.get('L9_DS', '100,150').split(',')]
DREF = float(os.environ.get('L9_DREF', '20000'))

if __name__ == '__main__':
    kind = sys.argv[1]
    for D in DS:
        for J in [int(a) for a in sys.argv[2:]]:
            nus = [x*P.mphys(D, J) for x in XS]
            P.get(kind, J, D, nus)
            P.get(kind, J, DREF, nus)
