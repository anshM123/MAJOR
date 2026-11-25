"""Shared setup for PARITY_ODD_DZ.md (ads4/code/oddz).  Imports code/kernels (pipeline, kernels_lib) and code/angproofs
read-only; KERNELS_CACHE is redirected to oddz/cache; every write is guarded to oddz/cache."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'cache')
os.environ['KERNELS_CACHE'] = CACHE
sys.path.insert(0, os.path.join(HERE, '..', 'kernels'))
import sympy as sp
import pipeline as pl          # read-only
import kernels_lib as kl       # read-only (cache redirected)
S = sp.Symbol('s'); Jsym = sp.Symbol('J')
ANGP = os.path.join(HERE, '..', 'angproofs')

def guard(path):
    p = os.path.abspath(path)
    assert p.startswith(CACHE + os.sep), 'refusing to write outside oddz/cache: ' + p
    return p

def write(name, lines):
    open(guard(os.path.join(CACHE, name)), 'w').write('\n'.join(lines) + '\n')

def symbolic_vectors():
    """ANGULAR_PROOFS Theorem 1 vectors (conserved + Bose symmetric for symbolic J), read-only from angproofs/cache.
       keys (m1, m3, k, n13, n14, n34, eps): monomial V1^m1 V3^m3 V4^(J-k) H13^n13 H14^n14 H34^n34 [eps_ij X_ij]."""
    import pickle
    return pickle.load(open(os.path.join(ANGP, 'cache', 'symbolic_structures.pkl'), 'rb'))

def structure_at_J(name, J):
    """(monos, coeffs) in pipeline monomial format (m1, m3, m4, n13, n14, n34, eps) at integer J; coeffs polynomials in s."""
    vec = symbolic_vectors()[name]
    monos, coeffs = [], []
    for (m1, m3, k, n13, n14, n34, eps), c in vec.items():
        c = sp.expand(sp.sympify(c).subs(Jsym, J))
        if c == 0 or J - k < 0:
            continue
        monos.append((m1, m3, J - k, n13, n14, n34, eps)); coeffs.append(c)
    return monos, coeffs

def moment_expr(expr, parity, mode):
    """first moments of the light-ray operator, identical to code/ladder/common.py (LADDER_ALLJ §1):
       'mid' : u-moment about the midpoint of x3, x4 (u = alpha/2): per monomial cf *= (a-b)/(2(a+b-2)) * (2 z.x) / 2 * 2;
       'orig': (u-moment about x3 = 0) - (midpoint moment) = -(z.x)/2 * I0.   a - b = 2s-5+i-j for both parities."""
    esum = 6 if parity == 'odd' else 5
    new = {}
    for (i, j, k, l, n1, n2, eps), cf in expr.items():
        if mode == 'mid':
            fac = (2*pl.s - 5 + i - j)*pl.C(1, 2*(esum - 2 - i - j))*2*pl.C(1, 2)
        elif mode == 'orig':
            fac = pl.C(-1, 2)
        key = (i, j, k, l + 1, n1, n2, eps)
        new[key] = new.get(key, pl.ZERO) + pl.redI(cf*fac)
    return {kk: v for kk, v in new.items() if v != 0}
