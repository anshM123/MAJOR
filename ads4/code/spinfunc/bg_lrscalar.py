"""
bg_lrscalar (BELOW_FOLD2.md §3.3 (c)): large-J approach of the resummed scalar expansion (bg_scalarJ, N = 12, R = 4) to the
leading-Regge image sum_j h_j P_{-1/2+i(nu-ij)}(cosh beta) (Stage-4 rank-1 scalar kernel sc1), along fixed bhat = (2J+1)/m:
  Gamma_LR(nu) = m^2 [log(R_bg/R_LR) - (same at nu_0)],  R = state/ref,  R_LR = T_st/T_ref.
It must converge to a finite function of (nu, bhat) as m grows (beyond-LR corrections are 1/m^2).
Usage: python3 bg_lrscalar.py -> cache/bg_lrscalar.txt
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
RQ = os.path.join(HERE, '..', 'resultq'); sys.path.insert(0, RQ)
import numpy as np
import bf_qan as BQ, bg_scalarJ as SJ
import stage4_hj as HJ4

if __name__ == '__main__':
    r = pickle.load(open(os.path.join(HERE, 'cache', 'bg_scalarJ_N12_R4.pkl'), 'rb'))
    nus = r['nus']
    F = HJ4.numfuncs('sc1'); R = {j: np.array([complex(F[j](complex(v))) for v in nus]) for j in F}
    lines = []
    for bh in (0.15, 0.3):
        lines.append(f'## bhat = (2J+1)/m = {bh}:  m^2[log(R_bg/R_LR)] - (nu_0 value), nu = ' + ' '.join(f'{v:g}' for v in nus))
        for D in (336.0, 1000.0, 3000.0, 10000.0):
            J = int(round(bh*(D - 2)/2)); J += J % 2
            m, eta = BQ.mass(D, J); _, etar = BQ.mass(2e4, J)
            if (2*J + 1)/D > 0.6:
                continue
            Rbg = SJ.ratio(r, D, J)
            Rlr = BQ.T(nus, eta, R)/BQ.T(nus, etar, R)
            G = m*m*np.log(Rbg/Rlr); G = G - G[0]
            lines.append(f'   Delta={D:7.0f} J={J:5d} m={m:8.1f}: ' + ' '.join(f'{x.real:+.2f}' for x in G))
    print('\n'.join(lines))
    open(os.path.join(HERE, 'cache', 'bg_lrscalar.txt'), 'w').write('\n'.join(lines) + '\n')
