"""nf4: heavy scan of Omega^{NF,SC} (NF02 uplift, no forward terms) -- same states, nodes, skeleton remainders and criteria as u4_scan.py.
Reuses the functional-independent skeleton caches cache/skel_u4_<Dg>_Df<Df>.pkl (computed on a miss).  Records -> cache/recs_nf4_<Dg>_Df<Df>.pkl.
Usage: OMEGAU_DF=3 python3 nf4_scan.py Dg [Dg ...]
"""
import sys, os, time, pickle
HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'omegau'); sys.path.insert(0, HERE); OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache')
import numpy as np
from nf_common import NF, TransformsNF, flat_parts_nf, action_nf, KS, SIGMA, state
from common import Transforms, gl01
from u4_scan import mfl_grid, jlist, NP
import remk

pn, pw = gl01(NP)

def run(Dg):
    Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0))
    F = NF(); TR = TransformsNF(F)
    pk = os.path.join(HERE, 'cache', f'skel_u4_{Dg}_{remk.TAG}.pkl')
    cache = pickle.load(open(pk, 'rb')) if os.path.exists(pk) else {}
    nmiss = 0
    W0 = Transforms({k: np.array([1.0]) for k in KS})
    recs = {h: [] for h in SIGMA}
    t0 = time.time()
    for mfl in mfl_grid():
        u = 1/mfl**2
        for h in SIGMA:
            js = jlist(h, mfl)
            if js.size == 0:
                continue
            Jm = int(js.max())
            pf, wf = gl01(max(400, 2*Jm + 100))
            Js, B2, B3, parts = flat_parts_nf(F, h, u, Jm, pf, wf)
            for J in js:
                st = state(Dg, Lam, mfl, int(J))
                if st is None:
                    continue
                key = (round(float(mfl), 6), int(J))
                if key not in cache:
                    _, W, P0 = W0.regge(Lam, st['beta'])
                    cache[key] = remk.skeleton_state(st, int(J), pn*np.sqrt(u), W, P0); nmiss += 1
                i = int(np.searchsorted(Js, J))
                r = action_nf(F, TR, h, st, int(J), mfl, Lam, pn, pw, cache[key], parts[i], B2[i], B3[i])
                r.update(mfl=float(mfl), J=int(J), m=st['m'])
                recs[h].append(r)
    if False:
        pickle.dump(cache, open(pk, 'wb'))
    pickle.dump(recs, open(os.path.join(OUT, f'recs_nf4_{Dg}_{remk.TAG}.pkl'), 'wb'))
    print(f"# nf4 Dg={Dg} {remk.TAG}: {sum(len(v) for v in recs.values())} states, {nmiss} new skeleton states, {time.time() - t0:.0f} s", flush=True)

if __name__ == '__main__':
    for Dg in [int(a) for a in sys.argv[1:]] or [100]:
        run(Dg)
