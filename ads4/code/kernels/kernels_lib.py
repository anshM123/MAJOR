"""
Driver utilities: conserved structures (cached), amplitudes per structure, norms, normalised amplitudes.
Conventions: see pipeline.py / structures.py.   s = (Delta + J)/2.
"""
import os, pickle, math
from fractions import Fraction as Fr
import sympy as sp
import mpmath as mp
import structures as st
import pipeline as pl

CACHE = os.environ.get('KERNELS_CACHE', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache'))
os.makedirs(CACHE, exist_ok=True)
S = sp.Symbol('s')

def _clean(vec):
    vec = [sp.together(sp.sympify(c)) for c in vec]
    den = sp.lcm_list([sp.denom(c) for c in vec])
    nums = [sp.cancel(c*den) for c in vec]
    g = sp.gcd_list([n for n in nums if n != 0])
    out = [sp.factor(sp.cancel(n/g)) for n in nums]
    # sign convention: first non-zero entry has positive leading coefficient
    for c in out:
        if c != 0:
            lead = sp.Poly(c, S).LC() if c.free_symbols else c
            if lead < 0:
                out = [-q for q in out]
            break
    return out

def conserved(J, parity, sval=None, seed=7):
    """monos, basis (list of coefficient lists; polynomials in s, or numbers if sval is given)"""
    tag = f"cons_J{J}_{parity}_{'sym' if sval is None else str(sval).replace('/', 'o')}_{seed}.pkl"
    fn = os.path.join(CACHE, tag)
    if os.path.exists(fn):
        return pickle.load(open(fn, 'rb'))
    monos, basis = st.conserved_structures(J, parity, s=(S if sval is None else sp.Rational(sval)), seed=seed)
    basis = [_clean(list(b)) for b in basis]
    if len(basis) > 1:
        M = sp.Matrix(basis)
        if sval is not None:
            M = M.rref()[0]
            basis = [_clean(list(M.row(i))) for i in range(M.rows)]
    pickle.dump((monos, basis), open(fn, 'wb'))
    return monos, basis

def sym_to_poly(c):
    """sympy polynomial in s (rational coefficients) -> flint poly"""
    c = sp.sympify(c)
    P = sp.Poly(c, S)
    out = pl.ZERO
    for (e,), co in P.terms():
        co = sp.Rational(co)
        out = out + pl.C(int(co.p), int(co.q))*pl.s**int(e)
    return out

def mono_amplitudes(J, parity, monos, source='plus', sig=None, n=(1, 0)):
    """list over monomials of flint polys (in t, d, s, I), and the common K"""
    tag = f"amp_J{J}_{parity}_{source}_{'sym' if sig is None else str(sig).replace('/', 'o')}_{str(n[0]).replace('/', 'o')}_{str(n[1]).replace('/', 'o')}_{len(monos)}.pkl"
    fn = os.path.join(CACHE, tag)
    if os.path.exists(fn):
        K, strs = pickle.load(open(fn, 'rb'))
        return [pl.CTX.from_dict(dct) for dct in strs], K
    kin = pl.Kin(source, n=n)
    exprs = [pl.scri_mono(m, kin) for m in monos]
    K = pl.Kbounds(exprs, parity)
    polys = [pl.amplitude(e, kin, parity, sig=sig, K=K)[0] for e in exprs]
    pickle.dump((K, [dict(p.to_dict()) for p in polys]), open(fn, 'wb'))
    return polys, K

def combine(polys, coeffs):
    tot = pl.ZERO
    for p, c in zip(polys, coeffs):
        if c != 0:
            tot = tot + (sym_to_poly(c) if not isinstance(c, type(pl.ZERO)) else c)*p
    return pl.redI(tot)

def poly_to_sym(p):
    """flint poly in (d, s, I) -> sympy"""
    d = sp.Symbol('d')
    out = 0
    for (et, ed, es, eI), c in p.to_dict().items():
        out += sp.Rational(int(c.p), int(c.q))*sp.Symbol("t")**int(et)*d**int(ed)*S**int(es)*sp.I**int(eI)
    return sp.expand(out)

def G_prefactor(parity, K, sv):
    """G of pipeline.amplitude at numerical s (mpmath)"""
    Ka, Kb, KB, Kg = K
    sv = mp.mpf(sv) if not isinstance(sv, mp.mpf) else sv
    h = mp.mpf(1)/2 if parity == 'odd' else mp.mpf(0)
    e13, e14, e34 = 5 - sv + h, sv + h, sv + h
    return (-4*mp.pi*1j*mp.exp(1j*mp.pi*e14)*mp.pi**2.5*mp.mpf(2)**(4 - 2*e34)
            * mp.rgamma(e14 + Ka)*mp.rgamma(e13 + Kb)*mp.rgamma(e34 + KB)*mp.rgamma(e34 - mp.mpf(1)/2 + Kg))

def norm_formula(Delta, J, m):
    """N_{O,m} for unit-normalised O (components O_m of O(zeta(a,b)) = sum a^{J+m} b^{J-m} O_m)"""
    return (mp.pi**2.5*mp.mpf(2)**(4 - 2*Delta - J)*mp.binomial(2*J, J + m)*mp.rf(Delta - J - 1, J - abs(m))
            * mp.rf(Delta - 1, abs(m))*mp.rgamma(Delta + J)*mp.rgamma(Delta - mp.mpf(1)/2))

N_T_UNIT = None
def norm_T_plus():
    """norm of T(e_+)|0> for C_T = 1:  16 N_{2,2}(Delta=3) = pi^2/3"""
    return 16*norm_formula(3, 2, 2)
