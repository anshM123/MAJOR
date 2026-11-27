"""md_probe.py -- boundary exponents on D' (pcyc coordinates) of the EE k=2 / scalar ratio Q2 = u'v' prefh/(m12^2 m34^2)
(= Q_k4 / R, R = ((v'-u')^2-1)/(4u'v')), to choose the product-quadrature weights for kind 'E2'.
Imports pcyc_dlib and rd_core read-only."""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'pcyc')); sys.path.insert(0, os.path.join(HERE, '..', 'rd'))
import numpy as np
import pcyc_dlib as P
import rd_core as C
np.seterr(all='ignore')


def Q2(g, Delta):
    Q, parts = C.Qratio(g, Delta, parts=True)
    up = g['m12']*g['m34']/(g['m13']*g['m24']); vp = g['m14']*g['m23']/(g['m13']*g['m24'])
    R = (parts['vmu']**2 - 1.0)/(4*up*vp)
    return Q, Q/R, R


def at(tx, rx, s, ty, ry, ph, al, Delta):
    o = P.point(*[np.atleast_1d(np.asarray(v, float)) for v in (tx, rx, s, ty, ry, ph, al)], Delta, want_Q=False)
    Q, Qb, R = Q2(o['g'], Delta)
    return o['L'], np.real(Q), np.real(Qb), np.real(R)


if __name__ == '__main__':
    Delta = float(sys.argv[1]) if len(sys.argv) > 1 else 10.0
    rng = np.random.default_rng(1)
    base = dict(tx=1.0, rx=0.4, s=0.3, ty=1.2, ry=0.55, ph=0.3, al=1.1)
    lim = {'tx': [1e-2, 1e-3, 1e-4], 'rx0': [1e-2, 1e-3, 1e-4], 'rx1': [1 - 1e-2, 1 - 1e-3, 1 - 1e-4],
           's1': [1 - 1e-2, 1 - 1e-3, 1 - 1e-4], 'ty': [1e-2, 1e-3, 1e-4], 'ry0': [1e-2, 1e-3, 1e-4], 'ry1': [1 - 1e-2, 1 - 1e-3, 1 - 1e-4],
           'ph': [np.pi/2 - 1e-2, np.pi/2 - 1e-3, np.pi/2 - 1e-4], 'al0': [1e-2, 1e-3, 1e-4], 'alpi': [np.pi - 1e-2, np.pi - 1e-3, np.pi - 1e-4]}
    key = {'tx': 'tx', 'rx0': 'rx', 'rx1': 'rx', 's1': 's', 'ty': 'ty', 'ry0': 'ry', 'ry1': 'ry', 'ph': 'ph', 'al0': 'al', 'alpi': 'al'}
    print(f'Delta = {Delta}: local exponents d log|.|/d log(dist) of scalar L, Q(k=4), Q2(k=2), R')
    for trial in range(3):
        b = {k: v*(1 + 0.3*rng.uniform(-1, 1)) for k, v in base.items()}
        print(f' base point {trial}: ' + ', '.join(f'{k}={v:.3f}' for k, v in b.items()))
        for nm, vals in lim.items():
            rows = []
            for v in vals:
                pt = dict(b); pt[key[nm]] = v
                L, Q, Qb, R = at(pt['tx'], pt['rx'], pt['s'], pt['ty'], pt['ry'], pt['ph'], pt['al'], Delta)
                dist = {'tx': v, 'rx0': v, 'rx1': 1 - v, 's1': 1 - v, 'ty': v, 'ry0': v, 'ry1': 1 - v, 'ph': np.pi/2 - v, 'al0': v, 'alpi': np.pi - v}[nm]
                rows.append((np.log(dist), L[0], np.log(abs(Q[0])), np.log(abs(Qb[0])), np.log(abs(R[0]))))
            rows = np.array(rows)
            sl = (rows[-1, 1:] - rows[-2, 1:])/(rows[-1, 0] - rows[-2, 0])
            print(f'   {nm:5s}: L {sl[0]:+6.2f}  Q {sl[1]:+6.2f}  Q2 {sl[2]:+6.2f}  R {sl[3]:+6.2f}')
    # large t (the B zero): growth rate in t
    print(' B-zero direction (t -> inf): d log|.|/dt')
    for nm in ('tx', 'ty'):
        rows = []
        for v in (20.0, 40.0):
            pt = dict(base); pt[nm] = v*Delta/10
            L, Q, Qb, R = at(pt['tx'], pt['rx'], pt['s'], pt['ty'], pt['ry'], pt['ph'], pt['al'], Delta)
            rows.append((pt[nm], L[0], np.log(abs(Q[0])), np.log(abs(Qb[0])), np.log(abs(R[0]))))
        rows = np.array(rows); sl = (rows[1, 1:] - rows[0, 1:])/(rows[1, 0] - rows[0, 0])
        print(f'   {nm}: L {sl[0]:+.4f}  Q {sl[1]:+.4f}  Q2 {sl[2]:+.4f}  R {sl[3]:+.4f}   (L+logQ2: {sl[0]+sl[2]:+.4f}; target -(1-5/Delta) = {-(1-5/Delta):+.4f})')
