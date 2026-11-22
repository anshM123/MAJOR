"""
rs2_parity (HRS_CLOSE.md §3, check B4'): the only monomial class that a non-real cycle could weigh differently from the real section
D is  sa odd with s, sf even  (the second, sa-odd cohomology class of the alpha factor; s-odd and sf-odd forms are exact).  In the
engine's EE pieces it is populated only by round-off (|c| <= 5e-13 of the largest coefficient; rs2_bookkeep).  This script bounds
its integrated weight: the worst case is a cycle carrying only one alpha-arc (upper arc 0 < alpha < pi), for which
  int_0^pi ca^r sa u^A v^B dalpha = sum_j C(r,j) (-1)^j 2^{A+j+B+1} B(A+j+1, B+1),
and the ratio  |sum over the class (upper arc)| / |the surviving image (full circle)|  is printed per order k at several nu.
Usage: python3 rs2_parity.py Delta J   -> cache/rs2_parity_D{Delta}_J{J}.txt
"""
import os, sys, time
os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.special import loggamma, comb
import bg_ser as S, bg_geom as GG, bg_ee as BE

C = os.path.join(HERE, 'cache')
NUS = np.array([0.02, 0.5, 1.02, 2.02, 3.02, 5.02, 8.02])


def lbeta(x, y):
    return loggamma(x) + loggamma(y) - loggamma(x + y)


def upper_arc(E, Cf, E0, E1, nus):
    """sum over monomials of int over (radial, s, phi full; alpha upper arc) for sa-odd monomials (E[:, 9] == 1)."""
    out = np.zeros(len(nus), complex)
    z = 1j*nus
    for e, c in zip(E, Cf):
        lg = 0
        for i in range(4):
            lg = lg + loggamma(1 + (e[i] + E0[i] + z*E1[i])/2)
        G = e[5] + E0[5] + z*E1[5]
        lg = lg + loggamma(0.5) + loggamma(G + 1) - loggamma(G + 1.5)
        Mv = e[6] + E0[6] + z*E1[6]
        lg = lg + loggamma(0.5) + loggamma((Mv + 1)/2) - loggamma(Mv/2 + 1)
        A = e[10] + E0[10] + z*E1[10]; B = e[11] + E0[11] + z*E1[11]
        acc = 0
        for j in range(int(e[8]) + 1):
            acc = acc + comb(int(e[8]), j)*(-1.0)**j*np.exp((A + j + B + 1)*np.log(2.0) + lbeta(A + j + 1, B + 1))
        out += c*np.exp(lg)*acc
    return out


def main(Dl, J):
    t0 = time.time()
    K = 4
    P = BE.ee_pieces(K + 1, Dl, J)
    Bf = GG.bfac(K)
    polys = [('vmu', P['vmu']), ('pref', P['pref'][0][0]), ('R', P['R'])]
    A = GG.assemble(P['logs'], polys, K)
    lines = [f'rs2_parity Delta={Dl} J={J} EE4 (0,0): class (sa odd, s even, sf even) on one alpha-arc vs the surviving image']
    for j, cj in A['corr'].items():
        Pk = S.lmul(cj, Bf.c, K)
        for k in range(K + 1):
            if not len(Pk[k]):
                continue
            E = Pk[k].E.astype(np.int64); Cf = Pk[k].C
            m = (E[:, 9] == 1) & (E[:, 4] == 0) & (E[:, 7] == 0)
            full = S.integrate_fast(E, Cf, A['E0'], A['E1'], NUS)
            bad = upper_arc(E[m], Cf[m], A['E0'], A['E1'], NUS) if m.any() else 0*NUS
            ratio = np.abs(bad)/np.maximum(np.abs(full), 1e-300)
            lines.append(f'  j={j} k={k}: {m.sum():6d} monomials; |class (one arc)| / |image| at nu = {list(NUS)}: '
                         + ' '.join(f'{x:.1e}' for x in ratio))
    lines.append(f'done {time.time() - t0:.0f}s')
    txt = '\n'.join(lines)
    print(txt, flush=True)
    open(os.path.join(C, 'rs2_parity_D%g_J%d.txt' % (Dl, J)), 'w').write(txt + '\n')


if __name__ == '__main__':
    main(float(sys.argv[1]), int(sys.argv[2]))
