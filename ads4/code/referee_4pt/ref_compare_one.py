"""Referee: compare independently regenerated mode caches with the authors' caches and recompute table entries."""
import os, sys, glob, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'fourpt')); sys.path.insert(0, os.path.join(HERE, 'e2e'))
from tmoments4 import TensorMoments4
import t5_analysis_l1 as t5
mine = sorted(glob.glob(os.path.join(HERE, "fourpt", "cache", "modes4", sys.argv[1])))
for f in mine:
    a = TensorMoments4.load(f)
    b = TensorMoments4.load(os.path.join(HERE, '..', 'fourpt', 'cache', 'modes4', os.path.basename(f)))
    dev = 0.0
    for ab in a.modes:
        for nu in (0.0, 3.0):
            ref = abs(a.moment(6, nu, ab, (1, -1, 1, -1))[0])
            for q in t5.QS:
                va, sa = a.moment(6, nu, ab, q); vb, sb = b.moment(6, nu, ab, q)
                v0, s0 = a.moment(6, nu, ab, (1, -1, 1, -1))
                dev = max(dev, abs(va*np.exp(sa - s0) - vb*np.exp(sb - s0))/abs(v0))
    print(f"{os.path.basename(f)}: max_(ab,nu,q) |mine - authors'| / |dominant mode| = {dev:.1e}")
    res = t5.run_case(a, ks=(6,), nus=(0.0, 1.0, 3.0))
    for nu in (0.0, 1.0, 3.0):
        line = f"  k=6 nu={nu}: "
        for name in ('v1v1', 'v2v2', 'v1v2', 'oo', 'ss', 'os'):
            r = res[(6, nu, name)]
            line += f"{name}: " + "/".join(t5.fmt(r['R'][c]) for c in t5.LEAD[name][:2]) + (f" S+={t5.fmt(r['Rp'])} S-={t5.fmt(r['Rm'])}" if name in ('v2v2', 'v1v2') else "") + "; "
        print(line)
