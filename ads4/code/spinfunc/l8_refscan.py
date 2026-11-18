"""
l8_refscan (M6 round 4): reference evaluations for removing the nu-only normalisation of C^_2 EXACTLY (doc 22.2):
for every physical point (Delta, x) of a state grid, the same nu = x m(Delta) is evaluated at heavy reference states Delta_ref
(default 2e4, 4e4), plus basis points (x = 0.2 ... 0.8) at Delta_ref for its helicity basis.  Also (optionally) new fixed-x state points.
Usage: python3 l8_refscan.py J  [env L8_DSTATE=400,800,1600,3200  L8_XS=0.1,...,0.8  L8_DREF=20000,40000  L8_NEWD=6400,10000]
Resumable (l8_points cache).
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import l8_points as P

def env(k, d):
    return [float(v) for v in os.environ.get(k, d).split(',') if v]

if __name__ == '__main__':
    J = int(sys.argv[1])
    DS = env('L8_DSTATE', '400,800,1600,3200'); XS = env('L8_XS', '0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8')
    DR = env('L8_DREF', '20000,40000'); ND = env('L8_NEWD', '')
    BX = env('L8_BX', '0.2,0.3,0.4,0.5,0.6,0.7,0.8')
    for D in ND:                                                    # new fixed-x state points
        P.get(J, D, [x*P.mphys(D, J) for x in XS])
    nus = sorted(set(round(x*P.mphys(D, J), 9) for D in DS + ND for x in XS))
    for Dr in DR:
        P.get(J, Dr, [x*P.mphys(Dr, J) for x in BX])                # basis points at the reference
        P.get(J, Dr, nus)
