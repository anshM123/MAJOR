"""LOWNU_SCOPE common tools.  Read-only imports from ../robust (rcommon).  Writes only to ./cache.
NF02 weights in exact rationals; dilation psi_{2,a}(p) = a^-1 psi2(p/a), psi_{3,a}(p) = a psi3(p/a) (support [0,a], p = nu/Lam)."""
import os, sys, json
from fractions import Fraction as Fr
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'cache')
sys.path.insert(0, os.path.join(HERE, '..', 'robust'))
import rcommon  # noqa: E402  (read-only)

def guard(path):
    p = os.path.abspath(path)
    assert p.startswith(CACHE + os.sep), 'refusing write outside lownu/cache: ' + p
    return p

def nf02_exact():
    d = json.load(open(os.path.join(HERE, '..', 'robust', 'functionals_noforward.json')))['NF02']
    p2 = [Fr(x) for x in d['psi2_exact']]; p3 = [Fr(x) for x in d['psi3_exact']]
    return p2, p3

def mom_exact(c, j):
    return sum(ci / (i + j + 1) for i, ci in enumerate(c))

def nf02_float():
    p2, p3 = nf02_exact()
    return np.array([float(x) for x in p2]), np.array([float(x) for x in p3])
