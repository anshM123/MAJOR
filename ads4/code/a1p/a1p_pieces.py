"""
a1p_pieces (A1PRIME.md section 5): new sigma-pieces for LIGHTSIDE_DIST's engine (T' = sigma x T_free, read-only lsd_core/lsd_mellin),
built exactly with sympy (Mellin amplitude M and d/d delta12 at fixed delta13, delta14 = 2 - delta12 - delta13).

  xsN = Q_s^N/(d12 - 1), xtN = Q_t^N/(d14 - 1)   s-/t-channel exchanges with a degree-2N residue (N = 1: EIN's spin-2 xs, xt).
         The pair XN = xsN + xtN is the 1<->3-symmetric object; its leading y^{2N-1} terms cancel, so J0(sigma) = 2N - 2,
         n_eff = J0 - 2 = 2N - 4 (LIGHTSIDE section 2).  The s-pole d12 = 1 (single-trace block) is taken term-wise.
  xu6 = Q_u^3/(d13 - 1)                         u-channel exchange, degree-6 residue: J0 = 6, n_eff = 4.
  c6  = s2^3                                    contact, Mellin degree 6: n_eff = 4 (the R^4 analogue of k = 4 in this family).
  top-Regge representatives (same top coefficient, nothing else at n >= k):
  r6  = D^6/8, r6u = D^6/(64 (d13 - 1)), rXN = c_N(d13) (D/2)^{2N-2}, D = d12 - d14,
        c_N(d13) = coefficient of y^{2N-2} of XN at large y = d12 (computed here).
"""
import os, sys
sys.dont_write_bytecode = True
import sympy as sp
import numpy as np

d12, d13 = sp.symbols('d12 d13')
d14 = 2 - d12 - d13
Qs = (d13 - d14)**2/4 + sp.Rational(1, 12)
Qt = (d13 - d12)**2/4 + sp.Rational(1, 12)
Qu = (d12 - d14)**2/4 + sp.Rational(1, 12)
s2 = d12**2 + d13**2 + d14**2
D = d12 - d14

EXPR = {'c6': s2**3, 'r6': D**6/8, 'xu6': Qu**3/(d13 - 1), 'r6u': D**6/(64*(d13 - 1)),
        'xu4': Qu**2/(d13 - 1), 'r4u': D**4/(16*(d13 - 1))}
STPOLE = {}          # piece -> residue polynomial Q^N at d12 = 1 (as a function of d13) for s-channel pieces
for N in (3, 4):
    EXPR[f'xs{N}'] = Qs**N/(d12 - 1)
    EXPR[f'xt{N}'] = Qt**N/(d14 - 1)
    STPOLE[f'xs{N}'] = (Qs**N).subs(d12, 1)
    y = sp.symbols('y')
    X = (EXPR[f'xs{N}'] + EXPR[f'xt{N}']).subs(d12, y)
    ser = sp.series(X, y, sp.oo, n=2*N + 2).removeO()          # large-y expansion
    lead = sp.Poly(sp.expand(sp.simplify(ser*y**0)), y) if False else None
    # coefficient extraction: X = sum_k a_k y^k + O(1/y); take top surviving power
    num, den = sp.fraction(sp.together(X))
    q, r = sp.div(sp.Poly(sp.expand(num), y), sp.Poly(sp.expand(den), y))
    co = q.all_coeffs()                                          # highest first
    deg = q.degree()
    assert sp.simplify(co[0]) != 0
    EXPR[f'rX{N}'] = sp.expand(co[0])*(D/2)**deg
    EXPR[f'_deg{N}'] = deg
    EXPR[f'_top{N}'] = sp.factor(co[0])
CONT = {'c6': (2/3, 2/3), 'r6': (2/3, 2/3), 'xu6': (0.35, 1.3), 'r6u': (0.35, 1.3), 'xu4': (0.35, 1.3), 'r4u': (0.35, 1.3),
        'xs3': (1.3, 0.35), 'xt3': (0.35, 0.35), 'xs4': (1.3, 0.35), 'xt4': (0.35, 0.35), 'rX3': (2/3, 2/3), 'rX4': (2/3, 2/3)}
PIECES = [k for k in CONT]
_f = {p: sp.lambdify((d12, d13), EXPR[p], 'numpy') for p in PIECES}
_df = {p: sp.lambdify((d12, d13), sp.diff(EXPR[p], d12), 'numpy') for p in PIECES}
_st = {p: sp.lambdify(d13, STPOLE[p], 'numpy') for p in STPOLE}


def M_amp(p, a, b):
    return _f[p](a, b) + 0*a + 0*b


def dM_amp(p, a, b):
    return _df[p](a, b) + 0*a + 0*b


def st_res(p, b):
    """single-trace pole residue polynomial at d12 = 1 (s-channel pieces only)."""
    return _st[p](b) if p in _st else None


if __name__ == '__main__':
    for N in (3, 4):
        print(f'X{N} = xs{N} + xt{N}: leading power y^{EXPR["_deg%d" % N]} (J0(sigma) = {EXPR["_deg%d" % N]}, n_eff = {EXPR["_deg%d" % N] - 2}); '
              f'top coefficient c_{N}(d13) = {EXPR["_top%d" % N]}')
    for p in ('xs3', 'xs4'):
        print(p, 'single-trace residue polynomial at d12 = 1:', sp.factor(STPOLE[p]))
    # self-check: dM vs finite difference
    a = 0.3 + 0.7j; b = 0.35 - 1.1j
    for p in PIECES:
        h = 1e-6
        fd = (M_amp(p, a + h, b) - M_amp(p, a - h, b))/(2*h)
        print(f'  {p}: |dM - FD|/|dM| = {abs(dM_amp(p, a, b) - fd)/abs(dM_amp(p, a, b)):.1e}')
