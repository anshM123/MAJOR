"""
po_k2cmp (PARITY_ODD_CLOSE.md sect. 2.2): independent cross-check of the parity-odd analytic-radial engine (po_ee_odd, EE k = 2 image)
against the project's DIRECT non-perturbative EE points at odd J (hc_points n3 path, cache/hc_pts_EE_J5_D400.pkl: the p8/w0 odd 1 x 1
integrand, state Delta = 400 and reference 2e4 on common random numbers).  Compared: m^2 [D(nu) - D(nu0)], D = log(state/reference).
Usage: python3 po_k2cmp.py  -> cache/po_k2cmp.txt (engine images cached in cache/po_eeq2_D*_J5.pkl)
"""
import os, sys, pickle
os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import po_ee_odd as PO
import bg_eeq as EQ

C = os.path.join(HERE, 'cache')


def img(D, J, nus, K=2, depth=3):
    fn = os.path.join(C, 'po_eeq2_D%g_J%d.pkl' % (D, J))
    if os.path.exists(fn):
        return pickle.load(open(fn, 'rb'))
    r = PO.images(D, J, K, depth, nus=np.asarray(nus), kinds=('EE2',))
    pickle.dump(r, open(fn, 'wb'))
    return r


if __name__ == '__main__':
    J, D = 5, 400.0
    db = pickle.load(open(os.path.join(C, 'hc_pts_EE_J5_D400.pkl'), 'rb'))
    nus = sorted(k[1] for k in db if k[0] == D and k[1] < 41)
    st = img(D, J, nus); rf = img(EQ.DREF, J, nus)
    m = np.sqrt((D - J - 2)*(D + J - 1))
    lines = [f'# po_k2cmp: J={J} Delta={D}, m={m:.2f}; engine EE k=2 (K=2 and K=1) vs direct n3 points; m^2[D(nu) - D(nu0)], nu0 = {nus[-1]:.3f}']
    Deng = {}
    for Kt in (1, 2):
        Is = EQ.value(st, 'EE2', 0, 0, Kt)[0]; Ir = EQ.value(rf, 'EE2', 0, 0, Kt)[0]
        Deng[Kt] = np.log(Is/Ir)
    Dd, ed = [], []
    for nu in nus:
        s, r = db[(D, nu)], db[(EQ.DREF, nu)]
        v = np.log(s['E'][0, 0]/r['E'][0, 0]) + s['logN'] - r['logN']
        Dd.append(v); ed.append(np.hypot(s['err'][0, 0]/abs(s['E'][0, 0]), r['err'][0, 0]/abs(r['E'][0, 0])))
    Dd = np.array(Dd); ed = np.array(ed)
    for i, nu in enumerate(nus):
        e2 = m*m*(Deng[2][i] - Deng[2][-1]); e1 = m*m*(Deng[1][i] - Deng[1][-1]); dd = m*m*(Dd[i] - Dd[-1])
        lines.append(f'  nu={nu:7.3f}: engine K2 {e2.real:+10.1f} (K1 {e1.real:+10.1f}) | direct {dd.real:+10.1f} +- {m*m*np.hypot(ed[i], ed[-1]):.0f}'
                     f' | engine - direct {(e2 - dd).real:+8.1f}  (Im: engine {e2.imag:+.1f}, direct {dd.imag:+.1f})')
    print('\n'.join(lines))
    open(os.path.join(C, 'po_k2cmp.txt'), 'w').write('\n'.join(lines) + '\n')
