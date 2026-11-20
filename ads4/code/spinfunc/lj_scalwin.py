"""
lj_scalwin (HLJ_ANALYTIC.md §4.1): the exact scalar skeleton (sm_exact caches) in the Bessel window nu < nu_B = c m/(2J+1) and beyond.
Gamma_s = m^2 (ex_2/E_2 - P_J) (c2 units, the convention of sm_offset's scalar control), minus the kappa image (sm_curvfix, C2 ++).
Prints max |Gamma_s - kappa| on nu < nu_B (c = 3, 5) and on the bulk, to compare with the window tolerances T_B of lj_contour.
Usage: python3 lj_scalwin.py
"""
import os, sys, pickle
os.environ['SM_CURV'] = 'none'; os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.special import j0
import sm_curvfix as CF
import ee1000_common as E
import remk

for Dg, tagpref in ((1000, None), (300, 'n0t')):
    db = pickle.load(open(os.path.join(HERE, 'cache', 'sm_exact_%d.pkl' % Dg), 'rb'))
    for key in sorted(db, key=lambda k: k[0]):
        if (tagpref is None and len(k := key) != 3) or (tagpref is not None and not (len(key) == 4 and key[3] == tagpref)):
            if not (tagpref is not None and key[0] == 24 and len(key) == 3):
                continue
        r = db[key]; st = r['st']; J = r['J']; m = st['m']
        x = np.asarray(r['x']); nu = np.asarray(r['nus'])
        b = 2*np.arcsinh((J + .5)/m)
        P = remk.Pders_vec(nu, np.cosh(b), nmax=0)[0]
        E2 = (2 - x**2)/(1 - x**2)**2
        Gs = m*m*(np.asarray(r['rho'][2]) + 2*(P - j0((2*J + 1)*x))/E2)
        y = Gs - CF.curv_image('C2', '++', J, m, x, E.Q.A.env, E.Q.A.SIG)
        out = []
        for c in (3, 5):
            nuB = c*m/(2*J + 1)
            sel = nu < nuB
            out.append('nu<%.1f (c=%d, %d nodes): max|G-k| %.2f' % (nuB, c, sel.sum(), np.abs(y[sel]).max() if sel.any() else np.nan))
        bulk = (nu >= 5*m/(2*J + 1)) & (x < 0.8)
        out.append('bulk nu>=nu_B(5), x<0.8: max|G-k| %.0f, rms %.0f' % (np.abs(y[bulk]).max(), np.sqrt(np.mean(y[bulk]**2))))
        print('Dg=%d J=%d %s: ' % (Dg, J, key[3:] if len(key) > 3 else '') + ' | '.join(out))
