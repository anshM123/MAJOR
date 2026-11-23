"""Edge structure, normalised to the envelope 2|c_l e^{l beta} F_b| instead of |Ghat_b| (removes oscillation nodes)."""
import os, sys, math
import mpmath as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cblocks
from A2_corner import c_l, F_b
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache', 'A2b_edge_envelope.txt')
out = []
for nu in (0.7, 3.1):
    Bs = {b: cblocks.recursion(nu, b, 200) for b in range(5)}
    for rho in (0.3, 0.85):
        row = []
        for beta in (2.0, 4.0, 6.0, 8.0, 10.0):
            r = rho * math.exp(-beta)
            Q = cblocks.scaled_legendre_vec([nu], r, math.cosh(beta), 202)
            m = 0
            for b in range(5):
                G = cblocks.Ghat_vec(Bs[b], Q, r)[0]
                mp.mp.dps = 30; l = mp.mpc(-0.5, nu)
                e = c_l(l) * mp.exp(l * beta) * F_b(rho, nu, b)
                m = max(m, abs(G - float(2 * mp.re(e))) / float(2 * abs(e)))
            row.append(m)
        s = f'nu={nu} rho={rho}: beta=2,4,6,8,10 -> max_b |Ghat_b - 2Re[edge]|/(2|edge|) = ' + '  '.join(f'{x:.1e}' for x in row)
        print(s); out.append(s)
open(OUT, 'w').write('\n'.join(out) + '\n')
