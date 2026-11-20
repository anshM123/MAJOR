"""
sm_foldfit (SMOOTH_BOUND.md §1.3): is the graviton correction in the fold window (nu = 10 ... 40) oscillating (prop. to D_h) or smooth?
Direct all-orders node values of EE_DG300_SMALLNU (z4 main run, cache/z4_chain_f4_besthold00.pkl: NL rows, Re Gamma_c = local
forward-normalised correction, estimator 'best'; Delta = 336 rows of J = 24, 28) minus the kappa curvature image (sm_curvfix), fitted by
c D_h + r and (c0 + c2 x^2) D_h + r on the direct nodes with NUMIN <= nu <= NUMAX (n3 contour, validated to |dc2| <~ 100-200 there).
Usage: python3 sm_foldfit.py [NUMIN NUMAX]  -> stdout
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ['SM_CURV'] = 'none'
import numpy as np
import sm_curvfix as CF
Q = CF.Q
lo = float(sys.argv[1]) if len(sys.argv) > 1 else 10.0
hi = float(sys.argv[2]) if len(sys.argv) > 2 else 40.0
d = pickle.load(open(os.path.join(HERE, 'cache', 'z4_chain_f4_besthold00.pkl'), 'rb'))
for (J, D), rec in sorted(d['nlrows'].items()):
    m = rec['m']; nus = np.asarray(rec['nus']); sel = np.asarray(rec['direct']) & (nus >= lo) & (nus <= hi)
    x = nus[sel]/m
    for h in ('++', '+-'):
        if (h, 'c') not in rec:
            continue
        y = np.asarray(rec[(h, 'c')]).real[sel] - CF.curv_image('C2', h, J, m, x, Q.A.env, Q.A.SIG)
        Dh = CF.Dsh(h, J, x)
        print(f"J={J} Delta={D:g} {h}: nodes nu = {np.round(nus[sel], 1).tolist()};  Re Gamma_c - kappa = {np.round(y).tolist()}")
        for basis, lab in (([Dh, np.ones_like(x)], 'c D + r'), ([Dh, x*x*Dh, np.ones_like(x)], '(c0 + c2 x^2) D + r')):
            A = np.array(basis).T
            c, *_ = np.linalg.lstsq(A, y, rcond=None); res = y - A @ c
            dof = max(len(y) - A.shape[1], 1); cov = np.linalg.inv(A.T @ A)*(res @ res)/dof
            print(f"     {lab:20s}: c = {c[0]:+7.1f} +- {np.sqrt(cov[0,0]):6.1f}   r = {c[-1]:+6.1f} +- {np.sqrt(cov[-1,-1]):5.1f}   "
                  f"(resid {np.sqrt(np.mean(res**2)):.0f} of {np.sqrt(np.mean(y**2)):.0f})")
