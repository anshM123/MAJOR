"""
bf_q (BELOW_FOLD.md): EE k = 4 (= k = 2 EE integrand x Lemma-B factor R) heavy images of low-J states on the real section D of
bf_core (no saddle), state / reference (Delta_ref = 2e4, same J) on common QMC points, plus the scalar k = 2 crystal on the same points
(validated against the exact Mellin skeleton).  Output: per-replica sums -> cache/bf_q_<tag>.pkl.

Usage:  python3 bf_q.py run Delta J m R [tag]       (m = log2 points per replica, R replicas; theta = 0)
        python3 bf_q.py scal Delta J                 (scalar state/ref vs exact, from the run cache)
"""
import os, sys, pickle, time
os.environ.setdefault('OMP_NUM_THREADS', '2')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import bf_core as B

DREF = 2e4
CPHI_G = 0.5                                           # phi-sampler exponent: the EE integrand is ~ cos^{1/2}(phi) at the y light cone
NUS = np.round(np.arange(0.02, 6.0, 0.04), 4)        # midpoints of the dnu = 0.04 de-shift grid, nu < 6
C = os.path.join(HERE, 'cache')


def fn(Dl, J, tag=''):
    return os.path.join(C, 'bf_q_D%g_J%d%s.pkl' % (Dl, J, tag))


def run(Dl, J, m, R, tag='', seed0=0):
    path = fn(Dl, J, tag)
    db = pickle.load(open(path, 'rb')) if os.path.exists(path) else dict(nus=NUS, reps={}, shifts={})
    from scipy.special import eval_legendre
    for r in range(R):
        seed = seed0 + r
        if seed in db['reps']:
            continue
        t0 = time.time()
        S = B.Sample(m, 1000 + seed, cphi=CPHI_G)
        rec = {}
        for D_ in (Dl, DREF):
            a = 1.0/D_
            P = B.grav_parts(S, D_, J, 'EE')
            for key in ('EE4', 'EE2'):
                w = P['pref'][0, 0]*(P['R'] if key == 'EE4' else 1.0)
                PE = dict(lA=P['lA'], lB=P['lB'], G=P['G'] + np.log(w + 0j))
                sh = db['shifts'].get((D_, key))
                v, sh = B.combine(PE, NUS, D_, S.logp, shift=sh); db['shifts'][(D_, key)] = sh
                rec[(D_, key)] = v
            tot = 0
            for n0 in (1, -1):
                Ps = B.scalar_parts(S, D_, 0, n0)
                if J:
                    Z, _, _ = S.raw(a, n0)
                    Ps['G'] = Ps['G'] + np.log(eval_legendre(J, B.t_ZW(Z).real) + 0j)
                sh = db['shifts'].get((D_, 'sc', n0))
                v, sh = B.combine(Ps, NUS, D_, S.logp, shift=sh); db['shifts'][(D_, 'sc', n0)] = sh
                tot = tot + v
            rec[(D_, 'sc')] = tot
        db = pickle.load(open(path, 'rb')) if os.path.exists(path) else db
        db['reps'][seed] = rec
        pickle.dump(db, open(path + '.tmp', 'wb')); os.replace(path + '.tmp', path)
        print(f'bf_q run D={Dl} J={J} seed={seed} 2^{m}: {time.time() - t0:.0f}s', flush=True)


def ratios(Dl, J, key, tag=''):
    """state/ref ratio per replica (complex array R x nu), using the common shifts."""
    db = pickle.load(open(fn(Dl, J, tag), 'rb'))
    rs = [db['reps'][s][(Dl, key)]/db['reps'][s][(DREF, key)] for s in sorted(db['reps'])]
    return db['nus'], np.array(rs)


def scal(Dl, J, tag=''):
    import bf_exact as BE, n3_exact as NX
    nus, rs = ratios(Dl, J, 'sc', tag)
    rec, rref = BE.get(Dl, J), BE.get(DREF, J)
    Rex = np.array([NX.X(nu, rec)[0]/NX.X(nu, rref)[0] for nu in nus])
    R = rs.mean(0); dR = rs.std(0, ddof=1)/np.sqrt(len(rs))
    L = np.log(R/Rex); c = np.mean(L[nus < 1.0])
    m = rec['m']
    return nus, m*m*(L - c), m*m*dR/np.abs(R)


if __name__ == '__main__':
    if sys.argv[1] == 'run':
        Dl, J, m, R = float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
        run(Dl, J, m, R, sys.argv[6] if len(sys.argv) > 6 else '')
    elif sys.argv[1] == 'scal':
        nus, e, d = scal(float(sys.argv[2]), int(sys.argv[3]))
        for v, x, y in zip(nus, e, d):
            if round(v/0.04 - 0.5) % 5 == 0:
                print(f'nu={v:5.2f}  m^2 eps = {x.real:+9.2f} {x.imag:+9.2f}i  +- {y:7.2f}')
