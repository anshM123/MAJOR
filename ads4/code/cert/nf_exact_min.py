"""Exact (rational) FF and N = |B2|+|B3| of NF02 at low J as polynomials in u (Lemma P), and their ratio on a fine u grid.
Read-only use of ads4/code/robust (functionals_noforward.json, rcert.mom)."""
import sys, os, json
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); RB = os.path.join(HERE, '..', 'robust'); sys.path.insert(0, RB)
D = json.load(open(os.path.join(RB, 'functionals_noforward.json')))['NF02']
psi2 = [Fr(c) for c in D['psi2_exact']]; psi3 = [Fr(c) for c in D['psi3_exact']]
def mom(a, j): return sum(c / (i + j + 1) for i, c in enumerate(a))
SIG = {'++': 1, '+-': -1}
def e_h(h, i): return i * (i - 1) if h == '++' else (i + 3) * (i + 4)
def parts(h, J):
    """B2(u), B3(u) as exact polynomials in u (Lemma P)."""
    JJ = J * (J + 1); kmax = J if h == '++' else J - 4
    B2 = {}; B3 = {}; prod = Fr(1)
    for k in range(kmax + 1):
        if k > 0: prod *= Fr(-(JJ - e_h(h, k)), k * k)
        B2[k] = B2.get(k, 0) + prod * 2 * mom(psi2, 2 * k)
        B2[k + 1] = B2.get(k + 1, 0) - prod * mom(psi2, 2 * k + 2)
        B3[k + 1] = B3.get(k + 1, 0) + prod * SIG[h] * mom(psi3, 2 * k)
    n = max(max(B2), max(B3)) + 1
    return [B2.get(k, Fr(0)) for k in range(n)], [B3.get(k, Fr(0)) for k in range(n)]
def ev(P, u):
    r = Fr(0)
    for c in reversed(P): r = r * u + c
    return r
if __name__ == '__main__':
    best = []
    for h in ['++', '+-']:
        for J in (range(0, 40, 2) if h == '++' else range(4, 40)):
            B2, B3 = parts(h, J)
            for i in range(0, 2001):
                u = Fr(i, 2000)
                b2, b3 = ev(B2, u), ev(B3, u)
                N = abs(b2) + abs(b3)
                if N == 0: continue
                best.append((float((b2 + b3) / N), h, J, float(u)))
    best.sort()
    for r in best[:12]: print('FF/N = %.8f  %s J=%d u=%.4f' % r)
    for h, J in [('++', 4), ('+-', 5), ('++', 2), ('+-', 4)]:
        B2, B3 = parts(h, J); b2, b3 = ev(B2, Fr(1)), ev(B3, Fr(1))
        print(h, J, 'u=1 exact FF/N =', (b2 + b3) / (abs(b2) + abs(b3)), '=', float((b2 + b3) / (abs(b2) + abs(b3))), ' B2=%.6g B3=%.6g' % (b2, b3))
