"""
b4_shadow (BATCH4_FOLD.md §1.1): the <TTO> shadows needed to turn the batch-4 direct points into chain records (hc_nodes.to_record) on
rows that have no K2 row at Delta_g = 2000: the seam group's even (2020, 20) and odd (2023, 23).  Exactly the functions the project's
row builders call (b5fold/b5_shadow.py, extfold/ef_shadow.py):
  even: spinfunc/l5_shadow_matrix.shadow_matrix(Delta, J, 0.5, ntrip=6, seed=7)  -> (b, logK, resid) under (Delta, J)
  odd : spinfunc/p2_shadow_odd.shadow_all(Delta, J) (defaults)                  -> (b, logK, resid) under (Delta, J, DEG)
Checks: the same calls reproduce gapnear's (2012, 12) and (2021, 21, 6) bit for bit; the reference shadows (2e4, J) of every batch-4 J
are looked up in the candidate lists hc_nodes uses (spinfunc/cache, gapnear rows).
Usage: python3 b4_shadow.py   -> cache/b4_shadow_k2.pkl, cache/b4_p2shadow.pkl, cache/b4_shadow.txt
"""
import os, sys, glob, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import b4_guard  # noqa: F401
os.environ['L1_DEG'] = '6'
import warnings
warnings.filterwarnings('ignore')
import numpy as np
SF = os.path.abspath(os.path.join(HERE, '..', 'spinfunc'))
GNR = os.path.abspath(os.path.join(HERE, '..', 'gapnear', 'cache', 'rows'))
B4C = os.path.join(HERE, 'cache')
sys.path.insert(0, SF)
cwd = os.getcwd(); os.chdir(SF)
import l5_shadow_matrix as SM
import p2_shadow_odd as P2
os.chdir(cwd)
lines = []


def pr(s=""):
    print(s, flush=True); lines.append(s)


def cands(odd):
    C = os.path.join(SF, 'cache')
    if odd:
        fs = ([os.path.join(C, 'p2_shadow.pkl')] + sorted(glob.glob(os.path.join(C, 'q2_p2shadow*.pkl'))) +
              sorted(glob.glob(os.path.join(C, 'w1_p2shadow*.pkl'))) + sorted(glob.glob(os.path.join(C, 'ee1000_p2shadow*.pkl'))) +
              sorted(glob.glob(os.path.join(C, 'hc_p2shadow*.pkl'))) + [os.path.join(GNR, 'gn_p2shadow.pkl')])
    else:
        fs = ([os.path.join(C, 'l9_shadow_k2.pkl')] + sorted(glob.glob(os.path.join(C, 'q2_shadow_k2*.pkl'))) +
              [os.path.join(C, 'u5_shadow_k2.pkl'), os.path.join(C, 'r1_shadow_k2.pkl')] +
              sorted(glob.glob(os.path.join(C, 'w1_shadow_k2*.pkl'))) + sorted(glob.glob(os.path.join(C, 'ee1000_shadow_k2*.pkl'))) +
              sorted(glob.glob(os.path.join(C, 'y1_shadow_k2*.pkl'))) + sorted(glob.glob(os.path.join(C, 'hc_shadow_k2*.pkl'))) +
              [os.path.join(GNR, 'gn_shadow_k2.pkl')])
    d = {}
    for f in fs:
        if os.path.exists(f):
            for k, v in pickle.load(open(f, 'rb')).items():
                d.setdefault(k, []).append(os.path.basename(f))
    return d


def comp(a, b):
    return max(float(np.max(np.abs(np.asarray(a[0]) - np.asarray(b[0])))), abs(complex(a[1]) - complex(b[1])))


if __name__ == '__main__':
    os.chdir(SF)
    ge = pickle.load(open(os.path.join(GNR, 'gn_shadow_k2.pkl'), 'rb')); go = pickle.load(open(os.path.join(GNR, 'gn_p2shadow.pkl'), 'rb'))
    t = time.time(); r = SM.shadow_matrix(2012.0, 12, 0.5, ntrip=6, seed=7)
    pr(f"check even (2012, 12) vs gapnear: max |d(b, logK)| = {comp((r['b'], r['logK']), ge[(2012.0, 12)]):.3e} [{time.time() - t:.1f}s]")
    t = time.time(); r = P2.shadow_all(2021.0, 21); k = (2021.0, 21, P2.T.DEG)
    pr(f"check odd {k} vs gapnear: max |d(b, logK)| = {comp((r['b'], r['logK']), go[k]):.3e} [{time.time() - t:.1f}s]")
    ev, od = {}, {}
    t = time.time(); r = SM.shadow_matrix(2020.0, 20, 0.5, ntrip=6, seed=7); ev[(2020.0, 20)] = (r['b'], r['logK'], r['resid'])
    pr(f"new even (2020, 20): resid {np.max(np.abs(r['resid'])):.2e} [{time.time() - t:.1f}s]")
    t = time.time(); r = P2.shadow_all(2023.0, 23); k = (2023.0, 23, P2.T.DEG); od[k] = (r['b'], r['logK'], r['resid'])
    pr(f"new odd {k}: b_OO = {complex(np.asarray(r['b'])[-1, -1]):.6e}, logK = {complex(r['logK']):.6f}, resid {r['resid']} [{time.time() - t:.1f}s]")
    pickle.dump(ev, open(os.path.join(B4C, 'b4_shadow_k2.pkl'), 'wb')); pickle.dump(od, open(os.path.join(B4C, 'b4_p2shadow.pkl'), 'wb'))
    ce, co = cands(False), cands(True)
    for J in (12, 16, 20, 21, 23, 27, 29):
        odd = J % 2 == 1
        for D in (2000.0 + J, 20000.0):
            key = (D, J, 6) if odd else (D, J)
            src = (co if odd else ce).get(key, []) + (['b4 (new)'] if key in (od if odd else ev) else [])
            pr(f"   shadow {key}: {src if src else 'MISSING'}")
    open(os.path.join(B4C, 'b4_shadow.txt'), 'w').write("\n".join(lines) + "\n")
