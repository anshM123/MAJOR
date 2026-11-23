"""
poserr/eng_orders.py (POS_ERRORS.md §2.3): per-order size of the analytic-radial 1/Delta engine (BELOW_FOLD2 / HRS_CLOSE) for the EE
k = 2 (EE2, Result D's C'_EE) and k = 4 (EE4, Result Q's C4) images at J = 0, read from the EXISTING caches (read-only):
  spinfunc/cache/bg_eeq_D{300,1000}_J0_K4.pkl  (state Delta = 300, 1000, K = 4; and the reference Delta = 2e4, K = 2)
and from new runs at larger Delta written only to poserr/cache/eng_D{Delta}_J0_K4.pkl (bg_eeq.images is called in-process; it writes
nothing).  Output: c2-unit size of each order, m^2 |T_k(nu)| / A0(nu), with A0 the local envelope (max over +-1 in nu) of |T_0|,
and the same after removing the best nu-flat multiple of T_0 from T_1 (the "common normalisation" part, HLJ_CLOSE §1.1).
Usage: python3 eng_orders.py [run Delta ...]   (run: compute new engine images; otherwise only analyse)
"""
import os, sys, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
SF = os.path.join(HERE, '..', 'spinfunc')
sys.path.insert(0, SF)
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
CA = os.path.join(HERE, 'cache')


def mphys(D, J=0):
    al = D - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


def load_all():
    recs = {}
    for D in (300, 1000):
        db = pickle.load(open(os.path.join(SF, 'cache', 'bg_eeq_D%d_J0_K4.pkl' % D), 'rb'))
        recs[float(D)] = db['st']
        recs[2e4] = db['ref']                     # identical reference in both files (K = 2)
    for f in sorted(os.listdir(CA)):
        if f.startswith('eng_D') and f.endswith('_J0_K4.pkl'):
            r = pickle.load(open(os.path.join(CA, f), 'rb'))
            recs[float(r['Delta'])] = r
    return recs


def pieces(rec, kind):
    import bg_eeq
    I, P = bg_eeq.value(rec, kind, 0, 0)
    return P                                        # (K+1) x nu, includes the common prefactor


def envelope(a, nus, w=1.0):
    out = np.empty_like(a)
    for i, v in enumerate(nus):
        out[i] = a[np.abs(nus - v) <= w].max()
    return out


def analyse(recs, out):
    for kind in ('EE2', 'EE4'):
        out.write('\n== %s (J = 0): c2-unit size of order k, m^2 |T_k| / env|T_0|  (max over nu in windows)\n' % kind)
        out.write('%9s %9s | %s\n' % ('Delta', 'm', '  '.join('k=%d[0,2) [2,10) [10,41)' % k for k in range(1, 5))))
        rows = []
        for D in sorted(recs):
            rec = recs[D]
            P = pieces(rec, kind); nus = rec['nus']; m = mphys(D)
            A0 = envelope(np.abs(P[0]), nus)
            cols = []
            for k in range(1, 5):
                if k < P.shape[0]:
                    v = m*m*np.abs(P[k])/A0
                    cols.append([v[nus < 2].max(), v[(nus >= 2) & (nus < 10)].max(), v[nus >= 10].max()])
                else:
                    cols.append([np.nan]*3)
            # common (nu-flat) part of order 1: c = median of Re(T1/T0) away from small |T0|
            msk = np.abs(P[0]) > 0.3*A0
            c = np.median((P[1][msk]/P[0][msk]).real)
            rem = m*m*np.abs(P[1] - c*P[0])/A0
            rows.append((D, m, cols, m*m*c, [rem[nus < 2].max(), rem[(nus >= 2) & (nus < 10)].max(), rem[nus >= 10].max()]))
            out.write('%9.0f %9.1f | %s\n' % (D, m, '  '.join('%9.3g %9.3g %9.3g' % tuple(cc) for cc in cols)))
        out.write('   common part of order 1 (m^2 * median Re T1/T0) and remainder m^2|T1 - c T0|/env|T0| ([0,2) [2,10) [10,41)):\n')
        for D, m, cols, cm, rem in rows:
            out.write('   Delta = %8.0f: m^2 c = %10.4g  (= %.3f m);  remainder %9.3g %9.3g %9.3g\n' % (D, cm, cm/m, *rem))
    return


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'run':
        import bg_eeq
        for D in sys.argv[2:]:
            D = float(D); p = os.path.join(CA, 'eng_D%g_J0_K4.pkl' % D)
            if os.path.exists(p):
                continue
            t = time.time()
            r = bg_eeq.images(D, 0, 4, 5, verbose=False)
            pickle.dump(r, open(p + '.tmp', 'wb')); os.replace(p + '.tmp', p)
            print('engine D = %g: %.0fs' % (D, time.time() - t), flush=True)
    recs = load_all()
    with open(os.path.join(CA, 'eng_orders.txt'), 'w') as f:
        analyse(recs, f)
    print(open(os.path.join(CA, 'eng_orders.txt')).read())


# ------------------------------------------------------------------------------------------------ spin J (bj engine caches, read-only)
BJ_FILES = ['bj_eeq_D20000_J28_d2.pkl', 'bj_eeq_D336_J28_d2.pkl', 'bj_eeq_D336_J32_d2.pkl', 'bj_eeq_D20000_J28_d3.pkl', 'bj_eeq_D20000_J32_d2.pkl', 'bj_eeq_D20000_J64_d2_p.pkl', 'bj_eeq_D20000_J100_d2_p.pkl',
            'bj_eeq_D336_J28_d3.pkl', 'bj_eeq_D336_J32_d3.pkl', 'bj_eeq_D1100_J100_d2_p.pkl']


def analyse_bj(out):
    import bj_ee
    out.write('\n== EE2, spin J (bj engine, conserved 2x2 basis): per-order c2 size, max over nu-windows [<2, 2-6, 6-12, 12-20] (nu >= 20 excluded: D-hat truncation)\n')
    out.write('   common = entry mean, nu-flat part removed (median Re T1/T0 of the entry mean); struct = max_entries |m^2 (T_k^{ac}/T_0^{ac} - T_k^{00}/T_0^{00})|\n')
    for fname in BJ_FILES:
        p = os.path.join(SF, 'cache', fname)
        if not os.path.exists(p):
            out.write('   %s: missing\n' % fname); continue
        r = pickle.load(open(p, 'rb'))
        key = sorted(r['res'].keys())[0]                    # (N, Kr): largest N first in sort? take max N
        key = max(r['res'].keys())
        rec = r['res'][key]
        D, J = r['Delta'], r['J']; m = mphys(D, J); nus = rec['nus']
        ents = sorted(k for k in rec['res'] if k[0] == 'EE2')
        P = {e[1:]: bj_ee.value(rec, 'EE2', e[1], e[2])[1] for e in ents}
        K = P[(0, 0)].shape[0] - 1
        win = [nus < 2, (nus >= 2) & (nus < 6), (nus >= 6) & (nus < 12), (nus >= 12) & (nus < 20)]   # nu >= 20: D-hat series unconverged (POS_ERRORS §2.3)
        out.write('   (Delta, J) = (%g, %d), m = %.1f, J^2/Delta = %.3f, (N, Kr) = %s, orders 0..%d\n' % (D, J, m, J*J/D, key, K))
        # absolute deviations over a common smooth envelope (no ratios: the entries have different zeros)
        A0 = envelope(np.max([np.abs(P[e][0]) for e in P], axis=0), nus, 1.5)
        cs = {}
        for k in range(1, K + 1):
            tot = sum(P[e][0] for e in P); totk = sum(P[e][k] for e in P)
            good = np.abs(tot) > 0.3*envelope(np.abs(tot), nus, 1.5)
            c = np.median((totk[good]/tot[good]).real)            # nu-flat common factor of order k (entry sum)
            dev = np.max([m*m*np.abs(P[e][k] - c*P[e][0])/A0 for e in P], axis=0)
            st = np.max([m*m*np.abs((P[e][k] - c*P[e][0]) - (P[(0, 0)][k] - c*P[(0, 0)][0])*P[e][0]/np.where(np.abs(P[(0, 0)][0]) > 0, P[(0, 0)][0], 1))/A0
                         for e in P], axis=0) if len(P) > 1 else 0*nus
            f = lambda v: ' '.join('%8.3g' % v[w].max() if w.any() else '     nan' for w in win)
            out.write('     k=%d: nu-flat common m^2 c = %10.4g (= %.3f m) | max_entry |T_k - c T_0|/A0 %s\n' % (k, m*m*c, m*c, f(dev)))

if __name__ == '__main__' and os.environ.get('ENG_BJ', '1') == '1':
    with open(os.path.join(CA, 'eng_orders_bj.txt'), 'w') as f:
        analyse_bj(f)
    print(open(os.path.join(CA, 'eng_orders_bj.txt')).read())
