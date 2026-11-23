"""
ef_merge (EXT_FOLD.md §2): build the merged direct-point sets (extfold/cache/pts/{old,base,all}/) for the near-gap states folded here,
from read-only sources, and report every key collision with its agreement.

Sources (later ones override earlier ones on a common (Delta, nu) key; every collision is compared and logged):
  spinfunc/cache (HC) -> neargap2/cache -> highspin/cache -> smallnu/cache (top level; n1chk/, n3chk/ are cross-check caches, not used)
  -> ext jobs (ads4/ext/unpacked/out/<tag>/<job>/hc_pts_*.pkl; job A: n1 bulk rows of J = 7, 23, 27, 31).
J = 29 (job B): the n3 small-nu records at nu = 10, 12.5, 15, 17.5 are POOLED over every DISTINCT RQMC stream.  z2_points.n3_block
draws replica r of a job with seed0 = Z2_SEED from the scrambled-Sobol stream seed0 + r, and its proposal (adapt_g) is deterministic, so
job s and job s + 1 share 7 of their 8 replicas (checked bit for bit here).  The 16 jobs Z2_SEED = 5001 ... 5016 therefore contain the
23 distinct streams 5001 ... 5023 (not 128), and smallnu's original record (Z2_SEED = 4242) the 8 streams 4242 ... 4249: 31 distinct
replicas in all.  They are pooled into a record with 8 SUPER-replicas (stream groups of 4, 4, 4, 4, 4, 4, 4, 3 in stream order, each
super-replica = the mean of its group, so the record's mean is a weighted mean of the 31 iid replicas and the chain's 8-replica
jackknife works unchanged); 'E'/'err' are recomputed from the super-replicas.  State and reference are pooled with the SAME grouping
(common random numbers preserved).
Usage: python3 ef_merge.py   -> cache/pts/{old,base,all}/hc_pts_*.pkl, cache/ef_merge.txt, cache/ef_streams_J29.pkl
"""
import os, sys, glob, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ef_guard  # noqa: F401
import numpy as np

CODE = os.path.join(HERE, '..')
EXT = os.path.join(HERE, '..', '..', 'ext', 'unpacked', 'out')
EFC = os.path.join(HERE, 'cache'); PTS = os.path.join(EFC, 'pts')
SRC = [('HC', os.path.join(CODE, 'spinfunc', 'cache')), ('NG2', os.path.join(CODE, 'neargap2', 'cache')),
       ('HS', os.path.join(CODE, 'highspin', 'cache')), ('SN', os.path.join(CODE, 'smallnu', 'cache'))]
STATES = [(7, 1007.0), (23, 1023.0), (27, 1027.0), (29, 1029.0), (31, 1031.0)]
JOBA = {7: 'A_bulk_J7', 23: 'A_bulk_J23', 27: 'A_bulk_J27', 31: 'A_bulk_J31'}
NU3 = [10.0, 12.5, 15.0, 17.5]
DREF = 20000.0
lines = []


def pr(s=""):
    print(s, flush=True); lines.append(s)


def load(f):
    return pickle.load(open(f, 'rb')) if os.path.exists(f) else {}


def same(a, b):
    """bitwise-equal replicas?"""
    return np.array_equal(np.asarray(a['reps']), np.asarray(b['reps']))


def merge_into(out, add, tag, kind, J):
    for k, v in add.items():
        if k in out:
            a = out[k]
            d = np.abs(np.asarray(a['E']) - np.asarray(v['E'])).max(); s = max(np.abs(np.asarray(a['err'])).max(), np.abs(np.asarray(v['err'])).max(), 1e-300)
            pr(f"   collision {kind} J={J} key ({k[0]:g}, {k[1]:.6f}) from {tag}: identical replicas {same(a, v)}; |dE| = {d:.3e} ({d/s:.2f} sigma)")
        out[k] = v


def ext_job_files(tag, kind, J, D):
    return sorted(glob.glob(os.path.join(EXT, tag, '*', 'hc_pts_%s_J%d_D%g.pkl' % (kind, J, D))))


def streams_J29(kind):
    """{nu: {D: {stream id: replica array}}} over smallnu's seed-4242 record and the 16 ext jobs; checks shared streams bit for bit."""
    D = 1029.0
    out = {nu: {D: {}, DREF: {}} for nu in NU3}
    nchk = nbad = 0
    srcs = [os.path.join(CODE, 'smallnu', 'cache', 'hc_pts_%s_J29_D1029.pkl' % kind)]
    srcs += sorted(glob.glob(os.path.join(EXT, 'B_n3_J29_s*', '*', 'hc_pts_%s_J29_D1029.pkl' % kind)))
    for f in srcs:
        P = load(f)
        for (DD, nu), r in P.items():
            if nu not in NU3 or 'seed0' not in r:
                continue
            s0 = int(r['seed0'])
            for i, rep in enumerate(np.asarray(r['reps'])):
                sid = s0 + i
                prev = out[nu][DD].get(sid)
                if prev is not None:
                    nchk += 1
                    if not np.array_equal(prev, rep):
                        nbad += 1
                else:
                    out[nu][DD][sid] = np.asarray(rep)
    return out, nchk, nbad


GROUPS = None


def pool_record(base, st, ids):
    """pooled record: 8 super-replicas (groups of the sorted stream ids), mean / err recomputed; metadata from base."""
    ng = 8
    parts = np.array_split(np.arange(len(ids)), ng)
    sup = np.array([np.mean([st[ids[i]] for i in p], axis=0) for p in parts])
    r = dict(base)
    r['reps'] = sup
    r['E'] = sup.mean(0)
    r['err'] = sup.std(0, ddof=1)/np.sqrt(ng)
    r['R'] = ng
    r['pooled_streams'] = list(ids)
    r['pooled_groups'] = [[ids[i] for i in p] for p in parts]
    r['method'] = str(base.get('method', '')) + ' POOLED (ext_fold, %d distinct streams)' % len(ids)
    return r


if __name__ == '__main__':
    for st in ('old', 'base', 'all'):
        os.makedirs(os.path.join(PTS, st), exist_ok=True)
    pr("# ef_merge: direct-point sets for the near-gap states (Delta_g = 1000)")
    streams = {}
    for J, D in STATES:
        for kind in ('EE', 'C3H'):
            fn = 'hc_pts_%s_J%d_D%g.pkl' % (kind, J, D)
            old = {}
            pr(f"# {fn}")
            for tag, d in SRC:
                P = load(os.path.join(d, fn))
                if P:
                    pr(f"   {tag}: {sum(1 for k in P if k[0] == D)} state points")
                merge_into(old, P, tag, kind, J)
            base = load(os.path.join(CODE, 'highspin', 'cache', fn))
            allp = dict(old)
            if J in JOBA:
                fs = ext_job_files(JOBA[J], kind, J, D)
                n0 = sum(1 for k in allp if k[0] == D)
                for f in fs:
                    merge_into(allp, load(f), 'ext:' + os.path.basename(os.path.dirname(f)), kind, J)
                pr(f"   ext {JOBA[J]}: {len(fs)} job files -> state points {n0} -> {sum(1 for k in allp if k[0] == D)}")
            if J == 29:
                S_, nchk, nbad = streams_J29(kind)
                streams[kind] = S_
                pr(f"   J=29 n3 streams ({kind}): shared-stream replicas compared bit for bit: {nchk}, differing: {nbad}")
                for nu in NU3:
                    ids = sorted(set(S_[nu][D]) & set(S_[nu][DREF]))
                    pr(f"     nu={nu:5.1f}: distinct streams {len(ids)} ({ids[0]}..{ids[7]} | {ids[8]}..{ids[-1]})")
                    for DD in (D, DREF):
                        allp[(DD, nu)] = pool_record(old[(DD, nu)], S_[nu][DD], ids)
            for st, P in (('old', old), ('base', base), ('all', allp)):
                if P:
                    pickle.dump(P, open(os.path.join(PTS, st, fn), 'wb'))
    pickle.dump(streams, open(os.path.join(EFC, 'ef_streams_J29.pkl'), 'wb'))
    open(os.path.join(EFC, 'ef_merge.txt'), 'w').write("\n".join(lines) + "\n")
