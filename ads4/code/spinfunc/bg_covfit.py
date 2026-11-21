"""
bg_covfit (BELOW_FOLD2.md §4.2): decomposition of the measured Gamma_K (bg_qan) into a covariant renormalisation of the HYB constant
plus a remainder:  Gamma_K(nu) = dgamma sum_j R_j (nu - ij)^2 + kappa N^_4 + rem(nu)  (least squares on nu < NUFIT, weights 1/sum|R_j|).
A covariant term projects to s_d = dgamma exactly in every gauge (Stage 4: sum_j int chi_d R_j (nu-ij)^2 = 1); the gauge spread of s_d
is the projection of rem.  Usage: python3 bg_covfit.py Delta J [NUFIT] -> appends cache/bg_covfit.txt
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
RQ = os.path.join(HERE, '..', 'resultq'); sys.path.insert(0, RQ)
import numpy as np
import bf_qan as BQ, bg_qan as GQ

if __name__ == '__main__':
    D, J = float(sys.argv[1]), int(sys.argv[2]); NUF = float(sys.argv[3]) if len(sys.argv) > 3 else 12.0
    d = pickle.load(open(os.path.join(HERE, 'cache', 'bg_qan_D%g_J%d.pkl' % (D, J)), 'rb'))
    K = max(d['rs']); r = d['rs'][K]
    nus = r['nus']; R = BQ.Rj(nus)
    cov = sum(R[j]*(nus - 1j*j)**2 for j in R); NH = r['NH']; w = 1/sum(np.abs(R[j]) for j in R)
    sel = nus < NUF
    A = np.array([cov.real[sel]*w[sel], NH.real[sel]*w[sel]]).T
    co, *_ = np.linalg.lstsq(A, r['GK'].real[sel]*w[sel], rcond=None)
    rem = r['GK'] - co[0]*cov - co[1]*NH
    W = BQ.weights()
    P = GQ.project(dict(nus=nus, GK=rem), W)
    Pc = GQ.project(dict(nus=nus, GK=cov), W)
    lines = [f'## Delta={D} J={J}: dgamma = {co[0]:+.4f}, kappa = {co[1]:+.3e} (fit nu < {NUF});  per-harmonic remainder |rem|/sum|R_j| max on nu<2.5 / 2.5-18: '
             f'{np.max(np.abs(rem*w)[nus < 2.5]):.3f} / {np.max(np.abs(rem*w)[(nus >= 2.5) & (nus < 18)]):.3f} c2']
    lines.append('   check int chi_d sum R_j (nu-ij)^2 per gauge: ' + ' '.join(f'{g}:{Pc[g][0]:+.4f}' for g in W))
    lines.append('   remainder projections s_d(rem): ' + ' '.join(f'{g}:{P[g][0]:+.4f}' for g in W))
    print('\n'.join(lines))
    open(os.path.join(HERE, 'cache', 'bg_covfit.txt'), 'a').write('\n'.join(lines) + '\n')
