"""
bh_k2scal (BELOW_FOLD3.md §2.1): the scalar calibration of the K2 Laplace machinery above the fold, Delta = 400, J = 0:
exact Mellin skeleton (n3_exact.X, bf_exact records) vs the analytic-radial expansion (bg_scalar, K = 6) vs the Laplace expansion about
the crystal saddle (l1_crystal.laplace, orders K0, K1, K2; L1_DEG = 6 as the w1 rows), as D(nu) = log(state/ref) (ref 2e4),
in c2 units m^2 [D(nu) - D(nu_0)], nu_0 = 39.85 (and 19.92).  Read-only imports only.
Usage: L1_DEG=6 python3 bh_k2scal.py  -> cache/bh_k2scal.txt, .pkl
"""
import os, sys, pickle
os.environ.setdefault('L1_DEG', '6')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import bg_scalar as BS, bf_exact as BE, n3_exact as NX, l1_crystal as Cr

D, J, DREF = 400.0, 0, 2e4
NUS = np.array([19.924984316, 25.0, 30.0, 39.849968632, 50.0, 60.0, 79.699937265, 100.0, 119.549905897, 159.399874529])

if __name__ == '__main__':
    out = dict(nus=NUS)
    rs, rr = BE.get(D, J), BE.get(DREF, J)
    m = rs['m']
    ex = []
    for nu in NUS:
        try:
            ex.append(np.log(NX.X(nu, rs)[0]/NX.X(nu, rr)[0]))
        except Exception as e:
            print('exact fails at', nu, e); ex.append(np.nan)
    out['ex'] = np.array(ex)
    r = BS.build(6, NUS, verbose=False)
    out['bg'] = {Kt: np.log(BS.ratio(r, D, J, Kt=Kt)[0]) for Kt in (2, 3, 4, 5, 6)}
    lap = {0: [], 1: [], 2: []}
    for nu in NUS:
        a, b = Cr.laplace(nu, D, J), Cr.laplace(nu, DREF, J)
        lap[0].append(a['logY0'] - b['logY0'])
        lap[1].append(a['logY0'] + np.log(1 + a['K1']) - b['logY0'] - np.log(1 + b['K1']))
        lap[2].append(a['logY0'] + np.log(1 + a['K1'] + a['K2']) - b['logY0'] - np.log(1 + b['K1'] + b['K2']))
        out.setdefault('K1m2', []).append(m*m*(a['K1'] - b['K1'])); out.setdefault('K2m2', []).append(m*m*(a['K2'] - b['K2']))
    out['lap'] = {k: np.array(v) for k, v in lap.items()}
    out['m'] = m
    pickle.dump(out, open(os.path.join(HERE, 'cache', 'bh_k2scal.pkl'), 'wb'))
    L = []
    for i0 in (3, 0):
        L.append(f'# scalar Delta=400 J=0, m={m:.2f}: m^2[D(nu) - D(nu0)], nu0 = {NUS[i0]:.2f};  columns exact | bg K=6 (bg-ex) | Laplace K0, K1, K2 (minus exact) | m^2 K1, m^2 K2 (state - ref)')
        f = lambda v: m*m*(v - v[i0])
        E = f(out['ex'])
        for i, nu in enumerate(NUS):
            L.append(f'  nu={nu:7.2f}: ex {E[i].real:+10.2f} | bg-ex {(f(out["bg"][6]) - E)[i].real:+8.3f} (K4 {(f(out["bg"][4]) - E)[i].real:+8.2f})'
                     f' | lap-ex K0 {(f(out["lap"][0]) - E)[i].real:+9.1f} K1 {(f(out["lap"][1]) - E)[i].real:+9.1f} K2 {(f(out["lap"][2]) - E)[i].real:+9.1f}'
                     f' | m2K1 {out["K1m2"][i].real:+9.1f} m2K2 {out["K2m2"][i].real:+9.1f}')
    print('\n'.join(L))
    open(os.path.join(HERE, 'cache', 'bh_k2scal.txt'), 'w').write('\n'.join(L) + '\n')
