"""deshift/kern.py -- shared kernels for DESHIFT_CLOSE.md (read-only use of stored data; nothing in resultq/ or spinfunc/ is imported).

Leading-Regge kernels, read as DATA (exact rationals, sympy strings):
  * Q, k = 4 (EE3):  R_j(nu), j = +-1, +-3       from resultq/functionals_q4_deshift.json ['EE3']['R_j']
                     (identical to resultq/cache/stage4_hj.pkl ['EE3']; checked in ab_solve.py)
  * Q, k = 5 (C5_pp): R_j(nu), j = 0, +-2        from resultq/functionals_q5_deshift.json ['R_j_k5_pp']
  * D, k = 2 (EE):   R_-(nu), R_+(nu)            from spinfunc/cache/u2_hclosed.pkl ['EE'] (sympy srepr / str)
Each rational is converted to exact Gaussian-rational coefficient lists and evaluated by Horner in python-flint ball arithmetic (acb),
or in numpy complex128.
"""
import os, json, pickle
import numpy as np
import sympy as sp
from flint import acb, arb, fmpq, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
ADS4 = os.path.abspath(os.path.join(HERE, '..', '..'))
CODE = os.path.join(ADS4, 'code')
CACHE = os.path.join(HERE, 'cache')
nu_s = sp.Symbol('nu')


def _q(x):
    x = sp.Rational(x)
    return fmpq(int(x.p), int(x.q))


class Rat:
    """exact rational function of nu with Gaussian-rational coefficients."""
    def __init__(self, expr):
        expr = sp.together(sp.sympify(expr, locals={'nu': nu_s, 'I': sp.I}))
        num, den = sp.fraction(expr)
        self.expr = expr
        self.num = [sp.nsimplify(c) for c in sp.Poly(sp.expand(num), nu_s).all_coeffs()]
        self.den = [sp.nsimplify(c) for c in sp.Poly(sp.expand(den), nu_s).all_coeffs()]
        self._prec = None
        self.np_num = np.array([complex(c) for c in self.num]); self.np_den = np.array([complex(c) for c in self.den])

    def _acb(self):
        if self._prec != ctx.prec:
            cv = lambda c: acb(arb(_q(sp.re(c))), arb(_q(sp.im(c))))
            self.a_num = [cv(c) for c in self.num]; self.a_den = [cv(c) for c in self.den]
            self._prec = ctx.prec
        return self.a_num, self.a_den

    def __call__(self, z):
        """acb evaluation."""
        an, ad = self._acb()
        n = acb(0)
        for c in an:
            n = n*z + c
        d = acb(0)
        for c in ad:
            d = d*z + c
        return n/d

    def np(self, z):
        z = np.asarray(z, complex)
        return np.polyval(self.np_num, z)/np.polyval(self.np_den, z)


def load_C4():
    d = json.load(open(os.path.join(CODE, 'resultq', 'functionals_q4_deshift.json')))['EE3']['R_j']
    return {int(j): Rat(v) for j, v in d.items()}


def load_C4_pkl():
    D = pickle.load(open(os.path.join(CODE, 'resultq', 'cache', 'stage4_hj.pkl'), 'rb'))['EE3']
    return {int(j): Rat(v) for j, v in D.items()}


def load_C5():
    d = json.load(open(os.path.join(CODE, 'resultq', 'functionals_q5_deshift.json')))['R_j_k5_pp']
    return {int(j): Rat(v) for j, v in d.items()}


def load_D():
    D = pickle.load(open(os.path.join(CODE, 'spinfunc', 'cache', 'u2_hclosed.pkl'), 'rb'))
    a, b = D['EE']
    ev = lambda s: sp.sympify(eval(s, {**sp.__dict__}) if isinstance(s, str) else s)
    return Rat(ev(a)), Rat(ev(b))


def load_xp():
    """Stage-7 extended-precision weights (functionals_q7_xprec.json): grid nu_k = (k+1/2)/25. Returns dict of mpmath-free decimal strings
    and float arrays."""
    d = json.load(open(os.path.join(CODE, 'resultq', 'functionals_q7_xprec.json')))
    out = {}
    for key in ('C4_value', 'C4_deriv'):
        e = d[key]; N = e['N']
        out[key] = dict(nu=(np.arange(N) + 0.5)/25.0, w=e['chi'], wf=np.array([float(x) for x in e['chi']]))
    e = d['C5']; N = e['N']
    out['C5'] = dict(nu=(np.arange(N) + 0.5)/25.0, w=list(zip(e['c_re'], e['c_im'])),
                     wf=np.array([float(a) for a in e['c_re']]) + 1j*np.array([float(b) for b in e['c_im']]))
    return out
