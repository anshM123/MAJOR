"""cl_mackorth: check that the n = 0 Mack polynomials (exactskel/mack.py, CHMRS conventions, d = 3) are orthogonal
with respect to the Mellin measure at the pole,  w(T) = Gamma(T/2)^2 Gamma((tau-T)/2)^2  on Re T = tau/2
(T = the B-hat label; T' = T + 2Df - tau).  If so, the spin decomposition of an s-channel residue R(T') is its
continuous-Hahn (orthogonal) expansion, and cl_box.py may build the spin basis by Stieltjes/Lanczos from w alone.
Read-only import of ../exactskel/mack.py.  Output: out/mackorth.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'exactskel'))
import mpmath as mp
import mack

mp.mp.dps = 30
out = []
for tau in (mp.mpf('2.5'), mp.mpf('3.7')):
    Js = list(range(0, 6))
    polys = []
    for J in Js:
        g = mack.Qpoly_ycoeffs(float(tau + J), J, [0], dps=40)[0]   # Qpoly^0 = sum g_i ymack^i, ymack = T - tau/2
        polys.append([mp.mpf(x) for x in g])
    def w(y):
        T = tau/2 + 1j*y
        return mp.gamma(T/2)**2*mp.gamma((tau - T)/2)**2
    def ev(g, y):
        ym = 1j*y
        return sum(c*ym**i for i, c in enumerate(g))
    G = [[mp.quad(lambda y: w(y)*ev(polys[i], y)*ev(polys[j], y), [-mp.inf, 0, mp.inf]) for j in range(len(Js))] for i in range(len(Js))]
    out.append('tau = %s : Gram matrix <Q_i Q_j> / sqrt(<Q_i Q_i><Q_j Q_j>) (should be the identity)' % tau)
    for i in range(len(Js)):
        out.append('  ' + ' '.join('%9.2e' % float(abs(G[i][j])/mp.sqrt(abs(G[i][i]*G[j][j]))) for j in range(len(Js))))
print('\n'.join(out))
open(os.path.join(HERE, 'out', 'mackorth.txt'), 'w').write('\n'.join(out) + '\n')
