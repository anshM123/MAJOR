"""LADDER_ALLJ: the ward_engine computation with SYMBOLIC J and m (and Delta), i.e. for all spins at once.
A generic local window F = sum_{delta=-2..2} f_delta A^{J-m-delta} B^{J+m+delta} a^4 phi_{2-m-delta}(z)  (source helicity +2; with
--minus: b^4 phi_{-2-m-delta}, helicity -2) is pushed through the same rest-frame K action as ward_engine.py; the coefficient of
A^{J-m}B^{J+m} is extracted.  Output: (1) -K_z F row m  vs  (W-18) row m;  (2) D^mu K_mu F row m  vs  the §8 recurrence row m.
Usage: python3 -B ward_symbolic.py [--minus]  -> cache/ward_symbolic[_minus].txt"""
import sys
import sympy as sp
from common import write
import ward_engine as we
from ward_engine import A, B, a, b, z0, z1, z2, D

J, m = sp.symbols('J m')
minus = '--minus' in sys.argv
f = {d: sp.Symbol(f'f{d}'.replace('-', 'm')) for d in range(-2, 3)}
h = -2 if minus else 2
def phi(k):
    return (z1 + sp.I*z2)**k*z0**(-2 - k)          # valid homogeneous extension for any (symbolic) k
we.phi = phi
tpol = b**4 if minus else a**4
F = sum(f[d]*A**(J - m - d)*B**(J + m + d)*tpol*phi(h - m - d) for d in range(-2, 3))

def row(expr):
    e = expr.subs({a: 0 if minus else 1, b: 1 if minus else 0}).subs({z0: 1, z1: 1, z2: 0})
    e = sp.expand(sp.powsimp(sp.expand(e*A**(m - J)*B**(-J - m)), force=True))
    e = sp.expand(sp.powsimp(e.subs(B, 1), force=True))
    return e

def row0(expr):
    e = row(expr)
    # remove any leftover A-dependent terms (other monomials)
    return sp.expand(sum(t for t in sp.Add.make_args(e) if not t.has(A)))

out = []
pr = lambda x: (print(x, flush=True), out.append(str(x)))
T = (D - 1)*(D - 2); C = J*(J + 1)
Kz = row0(we.Kb(F, [1, 1, 0], 1))
W18 = sp.expand((T + C - 2*m*m - 18)*f[0] - (D - m - 2)*(J + m + 1)*f[1] - (D + m - 2)*(J - m + 1)*f[-1])
pr(f'source helicity {h:+d}, symbolic J, m, Delta')
pr(f'  -K_z F (row m) = {sp.collect(sp.factor_terms(-Kz), list(f.values()))}')
pr(f'  (W-18) row m   = {sp.collect(W18, list(f.values()))}')
pr(f'  equal: {sp.expand(-Kz - W18) == 0}')
Kmu = [we.Kb(F, [1 if n == mu else 0 for n in range(3)], 1) for mu in range(3)]
cons = row0(sum(we.todorov(Kmu[mu])[mu] for mu in range(3)))
cons = sp.collect(sp.expand(cons), list(f.values()), sp.factor)
pr(f'  D^mu K_mu F (row m) = {cons}')
rec = sp.expand(m*(m - 4)*(T + C - 2*m*m)*f[0] - (m - 3)*(m - 4)*(D - m - 2)*(J + m + 1)*f[1] - m*(m - 1)*(D + m - 2)*(J - m + 1)*f[-1])
if not minus:
    pr(f'  equals -1/2 x [§8 recurrence row m]: {sp.expand(cons + rec/2) == 0}')
else:
    recm = rec.subs(m, -m).subs({f[1]: sp.Symbol('t'), f[-1]: f[1]}).subs(sp.Symbol('t'), f[-1])
    pr(f'  mirrored recurrence (m -> -m, f_{{m+-1}} -> f_{{m-+1}}): ratio = {sp.factor(sp.cancel(cons/recm)) if recm != 0 else None}')
write('ward_symbolic%s.txt' % ('_minus' if minus else ''), out)
