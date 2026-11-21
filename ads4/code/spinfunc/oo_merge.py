"""
oo_merge (HYPOTHESES_CLOSE2.md §3): pool independent thimble seeds of the same small-nu node into one z2-format record (state and
reference separately; seeds are independent, state and reference of one seed share their random numbers).  Replica i of the pooled
record = point-weighted mean of replica i of every seed (weights 2^mq), so the 8 pooled replicas stay independent and the z3 jackknife
applies unchanged.  logN / logNs must agree between seeds (same saddle normalisation; asserted to 1e-6).
Output: cache/z2<out>_pts_<kind>_J<J>.pkl (read through Z3_PREFER=x,<out>).
Usage: python3 oo_merge.py J nu out tag1,tag2,...   (tag '' = the main z2 cache; e.g. ",s2,oo3,oo4")
"""
import os, sys, pickle
import numpy as np
C = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache')

if __name__ == '__main__':
    J = int(sys.argv[1]); nu = round(float(sys.argv[2]), 9); out = sys.argv[3]; tags = sys.argv[4].split(',')
    for kind in ('EE', 'C3H'):
        res = {}
        for D in (336.0, 20000.0):
            recs = []
            for t in tags:
                f = os.path.join(C, 'z2%s_pts_%s_J%d.pkl' % (t, kind, J))
                d = pickle.load(open(f, 'rb')) if os.path.exists(f) else {}
                if (D, nu) in d:
                    recs.append((t, d[(D, nu)]))
            assert len(recs) >= 2, (kind, D, nu, [t for t, _ in recs])
            r0 = recs[0][1]
            for t, r in recs[1:]:
                assert abs(r['logN'] - r0['logN']) < 1e-6 and abs(r['logNs'] - r0['logNs']) < 1e-6, (t, r['logN'], r0['logN'])
            w = np.array([2.0**r['mq'] for _, r in recs]); w = w/w.sum()
            reps = sum(wi*r['reps'] for wi, (_, r) in zip(w, recs))
            repsS = sum(wi*r['repsS'] for wi, (_, r) in zip(w, recs))
            R = reps.shape[0]
            m = dict(r0)
            m.update(reps=reps, repsS=repsS, E=reps.mean(0), err=reps.std(0, ddof=1)/np.sqrt(R), Es=repsS.mean(),
                     errs=repsS.std(ddof=1)/np.sqrt(R), method='pooled thimble seeds ' + '+'.join(t or 'main' for t, _ in recs),
                     seeds=[t for t, _ in recs], weights=w.tolist())
            res[(D, nu)] = m
            rel = [float(np.abs(r['err'][0, 0]/r['E'][0, 0])) for _, r in recs]
            print(f"{kind} J={J} D={D:g} nu={nu}: seeds {[t or 'main' for t, _ in recs]} weights {np.round(w, 3).tolist()} rel err (0,0) "
                  f"{np.round(rel, 3).tolist()} -> pooled {abs(m['err'][0, 0]/m['E'][0, 0]):.3f}")
        f = os.path.join(C, 'z2%s_pts_%s_J%d.pkl' % (out, kind, J))
        db = pickle.load(open(f, 'rb')) if os.path.exists(f) else {}
        db.update(res)
        pickle.dump(db, open(f + '.tmp', 'wb')); os.replace(f + '.tmp', f)
