"""
n4_run (NEARGAP100_SMALLNU.md): crystal integrals at small nu (below the fold) as the sum of the A+ and A- Lefschetz thimbles of
n4_thimble (n+ = n- = 1), for a STATE and its REFERENCE (Delta_ref = 2e4) at the same nu on COMMON RANDOM NUMBERS (same proposal per
thimble, adapted on the state, same Sobol seeds): the state/reference ratio, which is all the graviton chain uses, is then much more
accurate than either value (the thimble Jacobian 'fingers' are common to both).
Stored per (Delta, nu):
  * raw: per thimble wh in (A+, A-): Es, repsS, logN (scalar, n1 normalisation incl. (2 pi)^{7/2}); per kind E, reps, logN (n1_run
    convention); prop, diag;
  * combined in n1_run / p8 format (E = E+ + E- exp(logN- - logN+), logN = logN+, reps replica-wise), so that n1_analysis.grid_direct /
    p8_analysis can use them unchanged (key (Delta, round(nu, 9))).
Files: cache/n4_pts_<kind>_J<J>.pkl (even, kinds C2, C3H), cache/n4_odd_pts_<kind>_J<J>.pkl (odd, P8.build kinds), cache/n4_scal.pkl
(scalar-only validation points, key (Delta, J, round(nu, 9))).
Usage: OMP_NUM_THREADS=1 python3 n4_run.py even J nu Delta_s [Delta_r]     /   odd J nu Delta_s   /   scal J nu D1,D2,...
"""
import os, sys, pickle, time, warnings
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
warnings.filterwarnings('ignore')
import numpy as np
import n4_thimble as N4
import l1_crystal as Cr

MQ = int(os.environ.get('N4_M', '9'))
RQ = int(os.environ.get('N4_R', '8'))
RHO1 = float(os.environ.get('N4_RHO1', '0.8'))
NUC = float(os.environ.get('N4_NUC', '10.0'))
SEED = 4242
DREF = 20000.0


def mphys(Delta, J):
    al = Delta - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


def _save(path, key, val):
    import fcntl
    with open(path + '.lock', 'w') as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        db = pickle.load(open(path, 'rb')) if os.path.exists(path) else {}
        db[key] = val
        tmp = path + '.tmp%d' % os.getpid()
        pickle.dump(db, open(tmp, 'wb')); os.replace(tmp, path)
        fcntl.flock(lk, fcntl.LOCK_UN)


_TR = {}


def tracks(Delta, J, nu, kinds, more=()):
    """saddle records at nu (A+, A-), from a track nu_c -> all of (nu, *more) descending; cached in memory per (Delta, J, kinds)."""
    key = (round(Delta, 9), J, tuple(kinds))
    if key in _TR and round(nu, 9) in _TR[key]:
        return _TR[key][round(nu, 9)], 0.0
    allnus = sorted(set([nu] + list(more)), reverse=True)
    nus = [NUC] + [v for v in allnus if v < NUC - 1e-9]
    t = time.time()
    recs = N4.track(Delta, J, nus, kinds=kinds)
    d = _TR.setdefault(key, {})
    for i, v in enumerate(nus):
        d[round(v, 9)] = {wh: recs[wh][i] for wh in ('A+', 'A-')}
    return d[round(nu, 9)], time.time() - t


def run_pair(J, nu, Ds, kinds, m=MQ, R=RQ, verbose=True, more=()):
    """thimble sums for the blocks Ds (first = state, adapts the proposals) at nu; returns {Delta: raw dict}."""
    out = {}; props = {}
    for D_ in Ds:
        t0 = time.time()
        rec, tt = tracks(D_, J, nu, kinds, more=(more if D_ == DREF else ()))
        raw = dict(nu=nu, Delta=D_, J=J, m=mphys(D_, J), kinds=kinds, mq=m, R=R, rho1=RHO1, nuc=NUC, th={})
        for wh in ('A+', 'A-'):
            o = N4.thimble(rec[wh], D_, J, kinds, rho1=RHO1, m=m, R=R, seed0=SEED, prop=props.get(wh))
            props.setdefault(wh, o['prop'])
            raw['th'][wh] = dict(S=o['S'], diag=o['diag'], prop=o['prop'], Phi0=rec[wh]['Phi0'], wc=rec[wh]['wc'],
                                 **{k: o[k] for k in kinds})
        raw['logbeta'] = Cr.log_beta(D_, J)
        out[D_] = raw
        if verbose:
            d = {wh: raw['th'][wh]['diag'] for wh in ('A+', 'A-')}
            print(f"J={J} Delta={D_:.4f} nu={nu:.4f}: scalar rel err A+ {abs(raw['th']['A+']['S']['errs']/raw['th']['A+']['S']['Es']):.3f} "
                  f"A- {abs(raw['th']['A-']['S']['errs']/raw['th']['A-']['S']['Es']):.3f}; nbad {d['A+']['nbad']}/{d['A-']['nbad']} "
                  f"nbad_g {d['A+']['nbad_g']}/{d['A-']['nbad_g']} of {d['A+']['n']}; track {tt:.0f}s, total {time.time() - t0:.0f}s", flush=True)
    return out


def combine_scalar(raw):
    """log Y = log(Y+ + Y-) (n1 normalisation, -log beta), replica-wise values."""
    a, b = raw['th']['A+']['S'], raw['th']['A-']['S']
    lp = a['logN'] - raw['logbeta']; lm = b['logN'] - raw['logbeta']
    reps = a['repsS'] + b['repsS']*np.exp(lm - lp)
    return lp, reps


def combine_kind(raw, k):
    a, b = raw['th']['A+'][k], raw['th']['A-'][k]
    f = np.exp(b['logN'] - a['logN'])
    reps = a['reps'] + b['reps']*f
    sa, sb = raw['th']['A+']['S'], raw['th']['A-']['S']
    repsS = sa['repsS'] + sb['repsS']*np.exp(sb['logN'] - sa['logN'])
    R = reps.shape[0]
    return dict(ljc=True, E=reps.mean(axis=0), err=reps.std(axis=0, ddof=1)/np.sqrt(R), reps=reps, logN=a['logN'],
                Es=repsS.mean(), errs=repsS.std(ddof=1)/np.sqrt(R), repsS=repsS, logNs=sa['logN'] - 3.5*np.log(2*np.pi),
                nu=raw['nu'], Delta=raw['Delta'], J=raw['J'], m=raw['m'], method='n4_thimble (A+ + A-)', mq=raw['mq'], R=R)


def fn(kind, J, odd):
    return os.path.join(HERE, 'cache', 'n4_%spts_%s_J%d.pkl' % ('odd_' if odd else '', kind, J))


if __name__ == '__main__':
    mode, J, nu = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
    if mode == 'scal':
        Ds = [float(v) for v in sys.argv[4].split(',')]
        res = run_pair(J, nu, Ds, ())
        for D_, raw in res.items():
            _save(os.path.join(HERE, 'cache', 'n4_scal.pkl'), (round(D_, 9), J, round(nu, 9)), raw)
    else:
        odd = mode == 'odd'
        Ds = [float(sys.argv[4]), float(sys.argv[5]) if len(sys.argv) > 5 else DREF]
        kinds = ('odd:C2', 'odd:C3H') if odd else ('C2', 'C3H')
        res = run_pair(J, nu, Ds, kinds)
        for D_, raw in res.items():
            _save(os.path.join(HERE, 'cache', 'n4_raw.pkl'), (mode, round(D_, 9), J, round(nu, 9)), raw)
            for k in kinds:
                kk = k.replace('odd:', '')
                _save(fn(kk, J, odd), (D_, round(nu, 9)), combine_kind(raw, k))
