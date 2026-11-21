"""
bg_validateJ (BELOW_FOLD2.md §3.2): scalar validation of the J/Delta-resummed expansion (bg_scalarJ, cache N12_R4 and N14_R2) against
the exact Mellin skeleton at large spin, convention of bg_validate (m^2 eps, const fitted on nu < 1.4).
Usage: python3 bg_validateJ.py -> cache/bg_validateJ.txt
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import bg_scalarJ as SJ, bf_exact as BE, n3_exact as NX

if __name__ == '__main__':
    r4 = pickle.load(open(os.path.join(HERE, 'cache', 'bg_scalarJ_N12_R4.pkl'), 'rb'))
    r2 = pickle.load(open(os.path.join(HERE, 'cache', 'bg_scalarJ_N14_R2.pkl'), 'rb'))
    nus = r4['nus']
    lines = ['# bg_validateJ: resummed (J/Delta) scalar expansion vs exact skeleton; m^2 eps per nu = ' + ' '.join(f'{v:g}' for v in nus)]
    for (D, J) in ((336., 24), (336., 28), (336., 32), (100., 24), (300., 24), (1000., 24), (100., 4), (300., 4)):
        rs, rr = BE.get(D, J), BE.get(2e4, J); m = rs['m']
        Rex = np.array([NX.X(nu, rs)[0]/NX.X(nu, rr)[0] for nu in nus])
        for lab, r, Nt, Rt in (('N12 R4', r4, 12, 4), ('N10 R4', r4, 10, 4), ('N14 R2', r2, 14, 2)):
            L = np.log(SJ.ratio(r, D, J, Nt=Nt, Rt=Rt)/Rex); L = m*m*(L - np.mean(L[nus < 1.4]))
            lines.append(f'  ({D:g},{J}) {lab}: max|.| nu<=4 / 4<nu<=10 (excl. 5-7) / 10-18: {np.abs(L.real[nus <= 4]).max():.2f} / '
                         f'{np.abs(L.real[(nus > 4) & (nus <= 10) & ((nus < 5) | (nus > 7))]).max():.2f} / {np.abs(L.real[nus > 10]).max():.2f} | '
                         + ' '.join(f'{x.real:+.2f}' for x in L))
    print('\n'.join(lines))
    open(os.path.join(HERE, 'cache', 'bg_validateJ.txt'), 'w').write('\n'.join(lines) + '\n')
