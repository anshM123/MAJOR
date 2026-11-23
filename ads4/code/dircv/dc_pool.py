"""
dc_pool (DIRCV.md §4): pool the per-wave (+ P_J twin + reference) records of several independent-seed runs at the same nodes into
8 equal "super-replicas" (common-seed pairing of state and reference preserved inside each group), so that the unchanged
pw_gamma/fw_gamma pipeline (8-replica jackknife, fw anchors with 8 replicas) runs on the pooled data.
Writes cache/pw_pts_J<J>_D<D>_pool.pkl and ..._poolP.pkl (tag '_pool'; use DC_LOCAL=_pool,_poolP in dc_pipe).
Usage: python3 dc_pool.py J Delta tag1,tag2,... [Ptags] [outtag]  (tags in spinfunc/cache or dircv/cache; Ptags: separate tag list
       for the P_J-twin file, e.g. with P_J-only runs of dc_ptwin; outtag default '_pool')
"""
import dc_guard, os, sys, pickle, copy
import numpy as np

NG = 8
INTERLEAVE = os.environ.get('DC_INTERLEAVE', '0') == '1'   # super-replica k = replicas i with i % 8 == k (grouping check)


def load(J, D, t, sfx=''):
    for base in (os.path.join(dc_guard.SPIN, 'cache'), dc_guard.CACHE):
        f = os.path.join(base, 'pw_pts_J%d_D%g%s%s.pkl' % (J, D, t, sfx))
        if os.path.exists(f):
            return pickle.load(open(f, 'rb'))
    return {}


def merge(recs):
    """recs: list of pw_points records of the same node/Delta (same modes); returns one record with NG super-replicas."""
    out = copy.deepcopy(recs[0])
    nm = len(out['per_mode'])
    Rtot = sum(len(r['per_mode'][0]['reps']) for r in recs)
    assert Rtot % NG == 0, Rtot
    gs = Rtot//NG
    for i in range(nm):
        l0 = recs[0]['per_mode'][i]['logN']
        reps = np.concatenate([np.asarray(r['per_mode'][i]['reps'])*np.exp(r['per_mode'][i]['logN'] - l0) for r in recs])
        sup = reps.reshape(gs, NG).mean(0) if INTERLEAVE else reps.reshape(NG, gs).mean(1)
        out['per_mode'][i].update(reps=sup, E=sup.mean(), err=sup.std(ddof=1)/np.sqrt(NG), logN=l0)
        for k in out['g']:
            g0 = recs[0]['g'][k][i]['logN']
            rg = np.concatenate([np.asarray(r['g'][k][i]['reps'])*np.exp(r['g'][k][i]['logN'] - g0) for r in recs])
            supg = rg.reshape((gs, NG) + rg.shape[1:]).mean(0) if INTERLEAVE else rg.reshape((NG, gs) + rg.shape[1:]).mean(1)
            out['g'][k][i].update(reps=supg, E=supg.mean(0), err=supg.std(0, ddof=1)/np.sqrt(NG), logN=g0)
    out['R'] = NG; out['pooled_from'] = Rtot
    return out


if __name__ == '__main__':
    J, D = int(sys.argv[1]), float(sys.argv[2]); tags = sys.argv[3].split(',')
    ptags = sys.argv[4].split(',') if len(sys.argv) > 4 else tags
    out = sys.argv[5] if len(sys.argv) > 5 else '_pool'
    for sfx in ('', 'P'):
        dbs = [d for d in (load(J, D, t, sfx) for t in (ptags if sfx == 'P' else tags)) if d]
        keys = set().union(*[set(d) for d in dbs])
        pooled = {}
        for key in sorted(keys):
            recs = [d[key] for d in dbs if key in d]
            pooled[key] = merge(recs)
            print(f"  {sfx or 'W'} {key}: pooled {[len(r['per_mode'][0]['reps']) for r in recs]} replicas -> {NG} x {pooled[key]['pooled_from']//NG}")
        pickle.dump(pooled, open(os.path.join(dc_guard.CACHE, 'pw_pts_J%d_D%g%s%s.pkl' % (J, D, out, sfx)), 'wb'))
