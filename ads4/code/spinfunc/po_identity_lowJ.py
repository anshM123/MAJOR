"""
po_identity_lowJ (PARITY_ODD_CLOSE.md sect. 1.1 (D)): po_identity's kappa test on the previous agent's K2 OO rows at J = 4-10
(cache/hc_p3_{EE,C3H}_J<J>.pkl, rows 400 / 1600 + references).  Read-only.  Usage: python3 po_identity_lowJ.py -> cache/po_identity_lowJ.txt
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import l8_analysis as L8

C = os.path.join(HERE, 'cache')
out = []
for J in (4, 6, 8, 10):
    for kind in ('EE', 'C3H'):
        fn = os.path.join(C, 'hc_p3_%s_J%d.pkl' % (kind, J))
        if not os.path.exists(fn):
            continue
        db = {(k[0], k[1]): v for k, v in pickle.load(open(fn, 'rb')).items() if k[2] == 1.0}
        for D in sorted(set(k[0] for k in db)):
            recs = {k: v for k, v in db.items() if k[0] == D}
            Ree = {k: dict(M0=v['M0'][:2, :2], M2=v['M2'][:2, :2], M4=v['M4'][:2, :2], logC=v['logC'], logK=v['logK'],
                           slog=v['slog'], sS=v['sS'], m=v['m']) for k, v in recs.items()}
            Ree_ee = Ree
            if kind == 'C3H':
                fe = os.path.join(C, 'hc_p3_EE_J%d.pkl' % J)
                de = {(k[0], k[1]): v for k, v in pickle.load(open(fe, 'rb')).items() if k[2] == 1.0 and k[0] == D}
                Ree_ee = {k: dict(M0=v['M0'][:2, :2], M2=v['M2'][:2, :2], M4=v['M4'][:2, :2], logC=v['logC'], logK=v['logK'],
                                  slog=v['slog'], sS=v['sS'], m=v['m']) for k, v in de.items()}
            B = L8.basis(J, D, Ree_ee) or L8.basis(J, D, Ree_ee, xs_min=0.0)
            if B is None:
                continue
            Ti = B['Ti']; T = B['T']
            ks = sorted(recs)
            M0 = recs[ks[0]]['M0']
            w = np.linalg.lstsq(M0[:2, :].T, M0[2, :], rcond=None)[0]
            kap = []
            for k in ks:
                r = recs[k]; tot = r['M0'] + r['M2'] + r['M4']
                kap.append(tot[2, 2]/(Ti @ tot[:2, :2] @ Ti.T)[0, 0])
            kap = np.array(kap); wt = w @ T
            out.append(f'{kind} J={J} Delta={D:g}: {len(ks)} nu; w.t+-/w.t++ = {abs(wt[1]/wt[0]):.1e}; kappa mean {np.mean(kap).real:.7g}, '
                       f'rel spread {np.abs(kap/np.mean(kap) - 1).max():.1e}')
print('\n'.join(out))
open(os.path.join(C, 'po_identity_lowJ.txt'), 'w').write('\n'.join(out) + '\n')
