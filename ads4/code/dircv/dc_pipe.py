"""
dc_pipe (DIRCV.md): Gamma - kappa per node at (J, Delta) in exactly the pw_gamma/fw_gamma conventions, WITH the jackknife replicas
(fw_gamma.analyse keeps only their spread).  Configurations: estimator EST in {dir, cv}, record set TAGS (pw_points tags: '_w' waves,
'_wP' P_J on the same CRN, '_none' = fw records only, or extra tags found in dircv/cache), KEEPFW (keep fw nodes at other nu).
Records from dircv/cache (new runs) take precedence over spinfunc/cache when DC_LOCAL=1 for a tag.
Output: cache/pipe_<J>_<D>_<est>_<tags>_<keepfw>.pkl  {nus, m, a, y[(kd,h)], yj[(kd,h)] (R x n), img}
Usage: python3 dc_pipe.py J:Delta est tags keepfw
"""
import dc_guard, os, sys, pickle, time
import numpy as np
import pw_gamma as PG, fw_gamma as FG, pw_points as PP

SPINC = PP.C
LOCALTAGS = set(t for t in os.environ.get('DC_LOCAL', '').split(',') if t)
_load0 = PP.load


def _load(J, D, sfx=''):
    if PP.TAG in LOCALTAGS:
        f = os.path.join(dc_guard.CACHE, 'pw_pts_J%d_D%g%s%s.pkl' % (J, D, PP.TAG, sfx))
        return pickle.load(open(f, 'rb')) if os.path.exists(f) else {}
    return _load0(J, D, sfx)


PP.load = _load


def run(J, D, est, tags, keepfw, pr=print):
    FG.EST = est
    PG.TAGS[:] = list(tags)
    PG.KEEPFW = keepfw
    PG.STATE[J] = D
    t0 = time.time()
    g0 = FG.row(J, D)
    R = len(next(iter(FG.fwpts('EE', J).values()))['reps'])
    gj = [FG.row(J, D, jk=i) for i in range(R)]
    m = g0['m']; xs = g0['xs']
    out = dict(J=J, D=D, m=m, nus=xs*m, xs=xs, a=g0['a'], est=est, tags=tuple(tags), keepfw=keepfw, R=R)
    for kd in ('C2', 'C3H'):
        for h in ('++', '+-'):
            if (kd, h) not in g0:
                continue
            img = FG.CF.curv_image(kd, h, J, m, xs, FG.Q.A.env, FG.Q.A.SIG)
            out[(kd, h)] = dict(G=np.asarray(g0[(kd, h)], float), img=img, y=np.asarray(g0[(kd, h)], float) - img,
                                yj=np.array([np.asarray(gg[(kd, h)], float) - img for gg in gj]),
                                aj=[gg['a'] for gg in gj])
    pr(f"pipe J={J} D={D:g} est={est} tags={tags} keepfw={keepfw}: nodes {np.round(xs*m, 2).tolist()} a={g0['a']} ({time.time()-t0:.0f}s)")
    return out


def fname(J, D, est, tags, keepfw):
    return os.path.join(dc_guard.CACHE, 'pipe_%d_%g_%s_%s_%d.pkl' % (J, D, est, '+'.join(tags) or 'none', int(keepfw)))


def get(J, D, est, tags, keepfw, redo=False):
    f = fname(J, D, est, tags, keepfw)
    if os.path.exists(f) and not redo:
        return pickle.load(open(f, 'rb'))
    o = run(J, D, est, tags, keepfw)
    pickle.dump(o, open(f, 'wb'))
    return o


if __name__ == '__main__':
    J, D = int(sys.argv[1].split(':')[0]), float(sys.argv[1].split(':')[1])
    est = sys.argv[2]; tags = [t for t in sys.argv[3].split(',') if t != 'none'] if sys.argv[3] != 'none' else ['_none']
    keepfw = sys.argv[4] == '1'
    get(J, D, est, tags, keepfw, redo=True)
