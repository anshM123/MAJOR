"""
b5_merge (BATCH5_FOLD.md §1): merged direct-point sets for the batch-5 near-gap states (Delta_g = 1000, own rows Delta = 1000 + J),
from read-only sources -- extfold/ef_merge.py with the source list extended by the batch-5 job pickles.

Sources (later override earlier on a common (Delta, nu) key; every collision compared and logged):
  spinfunc/cache (HC) -> neargap2/cache -> highspin/cache -> smallnu/cache (top level) -> [J = 23, 31: extfold/cache/pts/all, which
  already holds the job-A bulk rows and the batch-2 J = 31 GK15 edge] -> batch 5 (ads4/ext/unpacked/out/N5_J<J>_<grp>/<job>/hc_pts_*.pkl).
Sets written to b5fold/cache/pts/:
  all  : every point;   old : without batch 5;   base : highspin/cache only (sn_core's baseline normalisation, BASE[0]).
Also a per-job validity table (max relative replica error of state / reference, both kinds; > 0.2 = failed, hs_core.failed_x rule).
Usage: python3 b5_merge.py   -> cache/pts/{all,old,base}/hc_pts_*.pkl, cache/b5_merge.txt
"""
import os, sys, glob, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import b5_guard  # noqa: F401
import numpy as np

CODE = os.path.join(HERE, '..')
EXT = os.path.join(HERE, '..', '..', 'ext', 'unpacked', 'out')
B5C = os.path.join(HERE, 'cache'); PTS = os.path.join(B5C, 'pts')
SRC = [('HC', os.path.join(CODE, 'spinfunc', 'cache')), ('NG2', os.path.join(CODE, 'neargap2', 'cache')),
       ('HS', os.path.join(CODE, 'highspin', 'cache')), ('SN', os.path.join(CODE, 'smallnu', 'cache'))]
JS = (17, 18, 19, 20, 22, 23, 24, 25, 26, 28, 30, 31, 32)
DREF = 20000.0
lines = []


def pr(s=""):
    print(s, flush=True); lines.append(s)


def load(f):
    return pickle.load(open(f, 'rb')) if os.path.exists(f) else {}


def mphys(D, J):
    al = D - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


def relerr(r):
    return float(np.max(np.abs(r['err'])/np.maximum(np.abs(r['E']), 1e-300))) if np.all(np.isfinite(r['E'])) else np.inf


def merge_into(out, add, tag, kind, J):
    nc = 0
    for k, v in add.items():
        if k in out:
            a = out[k]; nc += 1
            d = np.abs(np.asarray(a['E']) - np.asarray(v['E'])).max()
            same = np.array_equal(np.asarray(a['reps']), np.asarray(v['reps']))
            pr(f"   collision {kind} J={J} ({k[0]:g}, {k[1]:.6f}) from {tag}: identical replicas {same}; |dE| = {d:.3e}")
        out[k] = v
    return nc


if __name__ == '__main__':
    for st in ('old', 'base', 'all'):
        os.makedirs(os.path.join(PTS, st), exist_ok=True)
    pr("# b5_merge: direct-point sets for the batch-5 states (Delta_g = 1000, rows Delta = 1000 + J)")
    for J in JS:
        D = 1000.0 + J; m = mphys(D, J)
        for kind in ('EE', 'C3H'):
            fn = 'hc_pts_%s_J%d_D%g.pkl' % (kind, J, D)
            old = {}
            for tag, d in SRC:
                merge_into(old, load(os.path.join(d, fn)), tag, kind, J)
            if J in (23, 31):
                merge_into(old, load(os.path.join(CODE, 'extfold', 'cache', 'pts', 'all', fn)), 'EF:all', kind, J)
            base = load(os.path.join(CODE, 'highspin', 'cache', fn))
            allp = dict(old)
            fs = sorted(glob.glob(os.path.join(EXT, 'N5_J%d_*' % J, '*', fn)))
            n0 = sum(1 for k in allp if k[0] == D)
            for f in fs:
                merge_into(allp, load(f), 'b5:' + os.path.basename(os.path.dirname(os.path.dirname(f))), kind, J)
            pr(f"# J={J} {kind}: old state points {n0}, batch-5 job files {len(fs)} -> {sum(1 for k in allp if k[0] == D)}")
            for st, P in (('old', old), ('base', base), ('all', allp)):
                if P:
                    pickle.dump(P, open(os.path.join(PTS, st, fn), 'wb'))
        # validity per batch-5 node
        PE = load(os.path.join(PTS, 'all', 'hc_pts_EE_J%d_D%g.pkl' % (J, D)))
        PC = load(os.path.join(PTS, 'all', 'hc_pts_C3H_J%d_D%g.pkl' % (J, D)))
        bad = []
        for grp in sorted(glob.glob(os.path.join(EXT, 'N5_J%d_*' % J))):
            g = os.path.basename(grp)
            for jd in sorted(glob.glob(os.path.join(grp, '*'))):
                nu = float(os.path.basename(jd).split('_nu')[1].split('_s')[0])
                ks = [k for k in PE if k[0] == D and abs(k[1] - nu) < 1e-6*max(nu, 1)]
                if not ks:
                    bad.append((g, round(nu/m, 4), 'missing')); continue
                nu_ = ks[0][1]
                e = max(relerr(P[(DD, nu_)]) if (DD, nu_) in P else np.inf for P in (PE, PC) for DD in (D, DREF))
                if e > 0.2:
                    bad.append((g, round(nu_/m, 4), f'{e:.2g}'))
        pr(f"   batch-5 nodes with replica rel. err > 0.2 (failed) or missing: {bad}")
    open(os.path.join(B5C, 'b5_merge.txt'), 'w').write("\n".join(lines) + "\n")
