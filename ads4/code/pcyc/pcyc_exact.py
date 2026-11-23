"""
pcyc_exact (PCYC_CLOSE.md section 1.4): END-TO-END check, at LOW Delta, that the physical scalar functional is the integral over D'.
The exact CHMRS/Mellin skeleton (exactskel.skelexact.ExactState via spinfunc/n4_exact.record and n3_exact.X, all imported read-only;
no project cache written: the records are computed here and stored in pcyc/cache) gives X(nu; Delta) = g C(nu+i) + g~ C(nu-i), the
literal CLS-D k = 2 functional on a scalar block up to nu-only and Delta-only factors (bg_validate's convention).  Compared with
    R_D'(nu) = I_S(Delta1, nu)/I_S(Delta2, nu)      (pcyc_dquad: exact integrand over D', the cut at the light cones of x5 included)
through  eps(nu) = log(R_D'/R_ex) - const  (const fitted on nu <= 2; the dropped Delta-only factors).  At Delta = 30 the region of D'
beyond raw radius 0.3 carries O(1) of the |integrand| weight (RD_CLOSE section 2.4), so this tests the chain itself, cut included, not
only its Regge expansion.
Usage: python3 pcyc_exact.py Delta1 Delta2 -> cache/pcyc_exact_D<1>_D<2>.txt
"""
import os, sys, pickle, time, glob
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
SPF = os.path.join(HERE, '..', 'spinfunc')
sys.path.insert(0, SPF); sys.path.insert(0, HERE)
import warnings
warnings.filterwarnings('ignore')
import numpy as np


def record(Delta):
    fn = os.path.join(HERE, 'cache', 'exact_D%g_J0.pkl' % Delta)
    if os.path.exists(fn):
        return pickle.load(open(fn, 'rb'))
    import n4_exact as N4E
    t = time.time()
    nus = np.concatenate([np.linspace(0.05, 2.0, 14), np.linspace(2.25, 24.0, 88)])
    r = N4E.record(float(Delta), 0, nus=nus)
    r['time'] = time.time() - t
    pickle.dump(r, open(fn, 'wb'))
    return r


def dq_best(Delta):
    import pcyc_compare as PC
    S = PC.dq(Delta, 'S')
    return S


def main(D1, D2):
    import n3_exact as NX
    r1, r2 = record(D1), record(D2)
    S1, S2 = dq_best(D1), dq_best(D2)
    L = [f'# pcyc_exact: D\' scalar ratio I_S({D1:g})/I_S({D2:g}) vs the exact CHMRS skeleton X({D1:g})/X({D2:g}); '
         f'records {r1.get("time", 0):.0f}s / {r2.get("time", 0):.0f}s']
    for k in (-1, -2):
        if len(S1) < -k or len(S2) < -k:
            continue
        g1, nus1, I1 = S1[k]; g2, nus2, I2 = S2[k]
        nus = np.array([v for v in nus1 if np.any(np.abs(nus2 - v) < 1e-9)])
        R = np.array([I1[np.argmin(np.abs(nus1 - v))]/I2[np.argmin(np.abs(nus2 - v))] for v in nus])
        Rex = np.array([NX.X(v, r1)[0]/NX.X(v, r2)[0] for v in nus])
        lg = np.log(R/Rex)
        c = np.mean(lg[nus <= 2.0])
        eps = lg - c
        L.append(f'  grids {g1} / {g2}:  const = {c.real:+.6f}{c.imag:+.6f}i')
        for v, e, a, b in zip(nus, eps, R, Rex):
            L.append(f'    nu = {v:5.2f}:  eps = {e.real:+.2e} {e.imag:+.2e}i   |R_D\'| = {abs(a):.6e}  arg R_D\' - arg R_ex = {np.angle(a/b):+.2e}')
    open(os.path.join(HERE, 'cache', 'pcyc_exact_D%g_D%g.txt' % (D1, D2)), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main(float(sys.argv[1]), float(sys.argv[2]))
