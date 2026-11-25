"""Shared helpers for ads4/code/ladder (LADDER_ALLJ.md).  Read-only use of code/kernels/{pipeline,structures,kernels_lib}.py;
KERNELS_CACHE is redirected to ./cache so nothing is written in code/kernels.  All writes guarded to ./cache."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'cache')
os.environ['KERNELS_CACHE'] = CACHE
sys.path.insert(0, os.path.join(HERE, '..', 'kernels'))
import sympy as sp
import pipeline as pl          # read-only
import kernels_lib as kl       # read-only (cache redirected)

S = sp.Symbol('s'); Dl = sp.Symbol('Delta'); a_, b_ = sp.symbols('a b')

def guard(path):
    p = os.path.abspath(path)
    assert p.startswith(CACHE + os.sep), 'refusing to write outside ladder/cache: ' + p
    return p

def write(name, lines):
    open(guard(os.path.join(CACHE, name)), 'w').write('\n'.join(lines) + '\n')

def moment_expr(expr, parity, mode):
    """mode 'mid': u-moment about the midpoint of x3,x4 (u = alpha/2):  cf *= (a-b)/(2(a+b-2)) (2 z.x)  -> l+1, cf*(a-b)/(a+b-2)/2*2
       mode 'orig': u-moment about x3 = 0 (T at the origin) minus the midpoint one, = (alpha_mid/2) I0 = -(z.x)/2 I0."""
    esum = 6 if parity == 'odd' else 5
    new = {}
    for (i, j, k, l, n1, n2, eps), cf in expr.items():
        if mode == 'mid':
            fac = (2*pl.s - 5 + i - j)*pl.C(1, 2*(esum - 2 - i - j))*2*pl.C(1, 2)   # (a-b)/(2(a+b-2)) * 2 (for 2z.x) * 1/2 (u=alpha/2)
        elif mode == 'orig':
            fac = pl.C(-1, 2)
        key = (i, j, k, l + 1, n1, n2, eps)
        new[key] = new.get(key, pl.ZERO) + pl.redI(cf*fac)
    return {kk: v for kk, v in new.items() if v != 0}
