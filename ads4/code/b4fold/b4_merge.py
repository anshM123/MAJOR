"""
b4_merge (BATCH4_FOLD.md §1.1): unpack the batch-4 external records (ext/results/ext_results_batch4.json, lossless; ext/unpack.py's
decoder, imported read-only) into b4fold/cache/unpacked/ -- ONLY the batch-4 job directories N4_* / N4S_* (the file also carries
earlier batches' records, which are not used here) -- and merge them into direct-point sets for the Delta_g = 2000 own rows
Delta = 2000 + J (J = 12, 16, 21, 27, 29; seam group J = 20, 23).
No earlier direct point exists on these rows; the two local pilot jobs of GAP_NEARGAP (gapnear/cache/pilot/N4_J29_{lo,edge}, the same
ext/point.py code path) are compared with the external records at the same keys (odd J: no basis ambiguity, so they must agree to the
RQMC determinism level) but are NOT merged (the external record is used).
Sets written to b4fold/cache/pts/all/ (b4_basis.py writes 'fix' = even J mapped to the local structure basis).
Also a per-node validity table (max relative replica error of state / reference, both kinds; > 0.2 = failed, hs_core.failed_x rule).
Usage: python3 b4_merge.py   -> cache/unpacked/out/N4*/..., cache/pts/all/hc_pts_*.pkl, cache/b4_merge.txt
"""
import os, sys, glob, json, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import b4_guard  # noqa: F401
import numpy as np
EXTD = os.path.abspath(os.path.join(HERE, '..', '..', 'ext'))
sys.path.insert(0, EXTD)
import unpack as UP                                    # read-only import (decoder only)
B4C = os.path.join(HERE, 'cache'); PTS = os.path.join(B4C, 'pts'); UNP = os.path.join(B4C, 'unpacked')
PILOT = os.path.join(HERE, '..', 'gapnear', 'cache', 'pilot')
JS = (12, 16, 20, 21, 23, 27, 29)
DREF = 20000.0
lines = []


def pr(s=""):
    print(s, flush=True); lines.append(s)


def mphys(D, J):
    al = D - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


def relerr(r):
    return float(np.max(np.abs(r['err'])/np.maximum(np.abs(r['E']), 1e-300))) if np.all(np.isfinite(r['E'])) else np.inf


if __name__ == '__main__':
    src = os.path.join(EXTD, 'results', 'ext_results_batch4.json')
    recs = json.load(open(src))
    pr(f"# b4_merge: {len(recs)} records in {os.path.basename(src)}")
    n = 0; logs = {}
    for r in recs:
        top = r['dir'].split('/')[1] if r['dir'].startswith('out/') else r['dir']
        if not (top.startswith('N4_') or top.startswith('N4S_')):
            continue
        d = os.path.join(UNP, r['dir']); os.makedirs(d, exist_ok=True)
        pickle.dump(UP.dec(r['data']), open(os.path.join(d, r['file']), 'wb')); n += 1
        if r.get('log'):
            logs[r['dir']] = r['log']
    pr(f"# unpacked {n} batch-4 records (N4_*, N4S_*) into cache/unpacked/out/; job logs {len(logs)}")
    json.dump(logs, open(os.path.join(B4C, 'b4_joblogs.json'), 'w'), indent=0)
    os.makedirs(os.path.join(PTS, 'all'), exist_ok=True)
    jobs = [l.split() for l in open(os.path.join(EXTD, 'jobs_batch4.txt')) if l.strip() and not l.startswith('#')]
    pr(f"# job list: {len(jobs)} job lines, {sum(len(j[5].split(',')) for j in jobs)} nodes")
    for J in JS:
        D = 2000.0 + J; m = mphys(D, J)
        for kind in ('EE', 'C3H'):
            fn = 'hc_pts_%s_J%d_D%g.pkl' % (kind, J, D)
            P = {}
            fs = sorted(glob.glob(os.path.join(UNP, 'out', 'N4*_J%d_*' % J, '*', fn)))
            for f in fs:
                for k, v in pickle.load(open(f, 'rb')).items():
                    if k in P:
                        pr(f"   collision {kind} J={J} {k}: identical {np.array_equal(np.asarray(P[k]['reps']), np.asarray(v['reps']))}")
                    P[k] = v
            pickle.dump(P, open(os.path.join(PTS, 'all', fn), 'wb'))
            pr(f"# J={J} {kind}: {len(fs)} job files -> {sum(1 for k in P if k[0] == D)} state points, {sum(1 for k in P if k[0] == DREF)} reference points")
        PE = pickle.load(open(os.path.join(PTS, 'all', 'hc_pts_EE_J%d_D%g.pkl' % (J, D)), 'rb'))
        PC = pickle.load(open(os.path.join(PTS, 'all', 'hc_pts_C3H_J%d_D%g.pkl' % (J, D)), 'rb'))
        # completeness against the job list, validity per node
        want = [(j[0], float(v)) for j in jobs if int(j[2]) == J for v in j[5].split(',')]
        bad, miss = [], []
        for g, nu in want:
            ks = [k for k in PE if k[0] == D and abs(k[1] - nu) < 1e-6*max(nu, 1)]
            if not ks:
                miss.append((g, round(nu/m, 4))); continue
            nu_ = ks[0][1]
            e = max(relerr(P[(DD, nu_)]) if (DD, nu_) in P else np.inf for P in (PE, PC) for DD in (D, DREF))
            if e > 0.2:
                bad.append((g, round(nu_/m, 4), f'{e:.2g}'))
        pr(f"   J={J}: job-list nodes {len(want)}, missing {miss}; failed (replica rel. err > 0.2): {bad}")
    # pilot comparison (odd J = 29: no basis ambiguity)
    pr("# pilot check (gapnear/cache/pilot, local, same code path) vs external, J = 29 row 2029:")
    for g in ('N4_J29_lo', 'N4_J29_edge'):
        for kind in ('EE', 'C3H'):
            for f in sorted(glob.glob(os.path.join(PILOT, g, '*', 'hc_pts_%s_J29_D2029.pkl' % kind))):
                PL = pickle.load(open(f, 'rb'))
                PX = pickle.load(open(os.path.join(PTS, 'all', 'hc_pts_%s_J29_D2029.pkl' % kind), 'rb'))
                for k, v in PL.items():
                    if k in PX:
                        a, b = np.asarray(v['E']), np.asarray(PX[k]['E'])
                        pr(f"   {g} {kind} {k[0]:g} nu={k[1]:.4f}: |E_loc - E_ext|/|E| = {np.abs(a - b).max()/np.abs(b).max():.2e}; "
                           f"rel. replica err {relerr(v):.2e} / {relerr(PX[k]):.2e}")
                    else:
                        pr(f"   {g} {kind} {k}: not in the external set")
    open(os.path.join(B4C, 'b4_merge.txt'), 'w').write("\n".join(lines) + "\n")
