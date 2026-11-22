"""
qf_jcheck (RESULT_Q_FWD_CLOSE.md sect. 3.1): is the fixed-order a = 1/Delta expansion adequate at J = 4 ... 13 near the gap, or is the J/Delta
resummation needed?  Scalar skeleton (where both exist): fixed order K = 4 and K = 6 (bg_scalar cache K6) vs the J/Delta-resummed
expansion N = 12, R = 4 (bg_scalarJ cache; validated against the exact skeleton to <= 0.04 c2 at (300, 24), (336, 24-32), BELOW_FOLD2 §2.2).
  m^2 [log(R_fixed/R_resum) - const(nu < 1.4)]  in c2, on the common nu nodes.  Read-only imports: bg_scalar, bg_scalarJ.
Usage: python3 qf_jcheck.py -> cache/qf_jcheck.txt
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import bg_scalar as SC, bg_scalarJ as SJ

if __name__ == '__main__':
    r4 = pickle.load(open(os.path.join(HERE, 'cache', 'bg_scalarJ_N12_R4.pkl'), 'rb'))
    nJ = r4['nus']
    r6 = SC.build(6, nus=nJ, verbose=False)            # fixed-order moments on the resummed cache's nodes (in memory)
    idx = list(range(len(nJ))); ok = np.ones(len(nJ), bool)
    nus = nJ
    lines = ['# qf_jcheck: scalar, fixed order (K = 4, 6) vs J/Delta-resummed (N12 R4); m^2 [log ratio - const(nu<1.4)] at nu = '
             + ' '.join(f'{v:g}' for v in nus)]
    states = [(304.03, 4), (306.03, 6), (307.03, 7), (308.03, 8), (309.03, 9), (311.03, 11), (313.03, 13), (441.85, 7), (441.87, 8),
              (441.89, 9), (1007.01, 7), (1009.01, 9), (1011.01, 11), (336., 24)]
    for D, J in states:
        al = D - J - 2.0; m = np.sqrt(al*(al + 2*J + 1))
        Rr = SJ.ratio(r4, D, J, Nt=12, Rt=4)[ok]
        out = []
        for K in (4, 6):
            Rf = SC.ratio(r6, D, J, Kt=K)[0][idx][ok]
            L = np.log(Rf/Rr); L = m*m*(L - np.mean(L[nus < 1.4]))
            out.append(L.real)
        lines.append(f'  ({D:g},{J:2d}) J^2/D={J*J/D:.3f}: K=4 max|.| nu<=4 / 4-10 / 10-18: {np.abs(out[0][nus <= 4]).max():.3f} / '
                     f'{np.abs(out[0][(nus > 4) & (nus <= 10)]).max():.3f} / {np.abs(out[0][nus > 10]).max():.3f};  K=6: '
                     f'{np.abs(out[1][nus <= 4]).max():.3f} / {np.abs(out[1][(nus > 4) & (nus <= 10)]).max():.3f} / {np.abs(out[1][nus > 10]).max():.3f}'
                     ' | K=4: ' + ' '.join(f'{x:+.2f}' for x in out[0]))
    # projection proxy: delta Gamma(nu) ~ T_st(nu) L(nu) (Gamma_exact ~ T_st m^2 eps), L = scalar K=4 fixed-minus-resummed (c2),
    # interpolated on the dnu = 0.04 grid (nu in [5, 7] excluded: zero of g + g~, ill-conditioned log), nu > 18: L(18) x (nu/18)^2;
    # delta s_d = 0.04 sum chi_d Re(T_st L)  [M: assumes the EE image has the scalar's J/Delta truncation error]
    import bf_qan as BQ, qf_qan as QQ
    W = BQ.weights()
    lines.append('# projection proxy delta s_d (H0 / SB) of the K = 4 fixed-order J/Delta error, delta Gamma = T_st L_scalar:')
    for D, J in states[:-1]:
        al = D - J - 2.0; m = np.sqrt(al*(al + 2*J + 1))
        Rr = SJ.ratio(r4, D, J, Nt=12, Rt=4)
        Rf = SC.ratio(r6, D, J, Kt=4)[0]
        L = np.log(Rf/Rr); L = (m*m*(L - np.mean(L[nus < 1.4]))).real
        keep = (nus < 5) | (nus > 7)
        B = QQ.base(D, J)
        g = B['nus']
        Lg = np.interp(g, nus[keep], L[keep]); Lg = np.where(g > 18, L[keep][-1]*(g/18.)**2, Lg)
        dG = B['Tst']*Lg
        out = []
        for gg in ('H0', 'SB'):
            nu, cv, cd = W[gg]; k = min(len(g), len(cd))
            out.append(0.04*np.sum(np.asarray(cd)[:k]*dG[:k].real))
        lines.append(f'  ({D:g},{J:2d}): delta s_d = {out[0]:+.4f} (H0) / {out[1]:+.4f} (SB)')
    print('\n'.join(lines))
    open(os.path.join(HERE, 'cache', 'qf_jcheck.txt'), 'w').write('\n'.join(lines) + '\n')
