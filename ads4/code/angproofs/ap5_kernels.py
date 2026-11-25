"""Step 5 (proof): exact simplification of the symbolic-J extreme amplitudes of ap4 to closed form, and the
absolute kernels.  J = 2L (even structures) or 2L+1 (odd J); all Gamma functions are reduced to a few base Gammas by
integer shifts and the Legendre duplication formula, so that the ratio to the claimed closed form is a rational function
of (s, L) (shown to be 1) and the 'absent-sector' extreme amplitude is shown to vanish identically."""
import pickle, sys
from common import *
L = sp.Symbol('L', integer=True, nonnegative=True)
GS = {}

def gsym(base):
    if base not in GS:
        GS[base] = sp.Symbol('G[%s]' % sp.sstr(base))
    return GS[base]

def norm_gamma(arg):
    arg = sp.expand(arg)
    c = arg.as_coeff_Add()[0] if arg.as_coeff_Add()[0].is_number else 0
    lin = sp.expand(arg - c)
    if lin == 0:
        return sp.gamma(arg)
    if sp.expand(lin - 2*S) == 0:                     # duplication
        return (sp.Integer(2)**(arg - 1)/sp.sqrt(sp.pi)*norm_gamma(arg/2)*norm_gamma(arg/2 + sp.Rational(1, 2)))
    fl = sp.floor(c)
    base = lin + (c - fl)
    m = int(fl)
    G = gsym(base)
    if m >= 0:
        return G*sp.Mul(*[base + q for q in range(m)])
    return G/sp.Mul(*[base + m + q for q in range(-m)])

def reduce_expr(e):
    e = e.replace(sp.rf, lambda a, b: sp.gamma(a + b)/sp.gamma(a)).replace(sp.ff, lambda a, b: sp.gamma(a + 1)/sp.gamma(a - b + 1))
    e = e.replace(sp.factorial, lambda a: sp.gamma(a + 1)).replace(sp.binomial, lambda a, b: sp.gamma(a + 1)/(sp.gamma(b + 1)*sp.gamma(a - b + 1)))
    e = e.replace(sp.gamma, norm_gamma)
    return e

def to_rational(e):
    """powers with symbolic exponents -> combine; returns simplified expression"""
    e = sp.expand_power_exp(sp.expand_power_base(e, force=True))
    e = sp.powsimp(e, force=True, combine='exp')
    e = e.rewrite(sp.exp)
    e = sp.powsimp(sp.expand_power_exp(e), force=True)
    return sp.simplify(e)

def claimed(name, which, Jsub):
    """M/U claimed, U = -4 pi i e^{i pi e14} pi^{5/2};  KERNELS_D3 §5.3 closed forms x LADDER_ALLJ ladder factor"""
    h = S - Jsub; Delta = 2*S - Jsub
    U = -4*sp.pi*sp.I*sp.pi**sp.Rational(5, 2)
    if name == 'pp' and which == 'minusJ':    # M_0 (4)_J/(Delta-1-J)_J, M_0 = 3i(-1)^{J/2} e^{i pi s} pi^3 2^{3-J} (Delta-J-1)_J/(G(s)G(Delta+3)G(1-h))
        M = 3*sp.I*(-1)**L*sp.exp(sp.I*sp.pi*S)*sp.pi**3*sp.Integer(2)**(3 - Jsub)*sp.rf(4, Jsub)/(sp.gamma(S)*sp.gamma(Delta + 3)*sp.gamma(1 - h))
        return M/(U*sp.exp(sp.I*sp.pi*S))
    if name == 'odd' and which == 'minusJ':
        M = 3*(-1)**L*sp.exp(sp.I*sp.pi*S)*sp.pi**3*sp.Integer(2)**(4 - Jsub)*sp.rf(4, Jsub)/(sp.gamma(S + sp.Rational(1, 2))*sp.gamma(Delta + 3)*sp.gamma(sp.Rational(3, 2) - h))
        return M/(U*sp.exp(sp.I*sp.pi*(S + sp.Rational(1, 2))))
    if name == 'pm' and which == 'plusJ':     # M_4 (4)_{J-4}/(Delta-1-J)_{J-4}
        M = -3*sp.I*(-1)**L*sp.exp(sp.I*sp.pi*S)*sp.pi**3*sp.Integer(2)**(2 - Jsub)*sp.rf(4, Jsub - 4)/(sp.gamma(S)*sp.gamma(Delta - 1)*sp.gamma(1 - h))
        return M/(U*sp.exp(sp.I*sp.pi*S))
    if name == 'oddJ' and which == 'plusJ':   # (-1)^{(J-5)/2} = (-1)^{L-2} = (-1)^L
        M = -3*sp.I*(-1)**L*sp.exp(sp.I*sp.pi*S)*sp.pi**3*sp.Integer(2)**(4 - Jsub)*sp.rf(4, Jsub - 4)/(sp.gamma(S + sp.Rational(1, 2))*sp.gamma(Delta - 1)*sp.gamma(sp.Rational(3, 2) - h))
        return M/(U*sp.exp(sp.I*sp.pi*(S + sp.Rational(1, 2))))
    return 0

TW, SG = sp.symbols('TW SG')     # TW = 2^L, SG = (-1)^L

def norm_powers(e):
    """replace every power / exponential with L in the exponent by monomials in TW = 2^L, SG = (-1)^L"""
    def fix_pow(b, ex):
        ex = sp.expand(ex)
        if L not in ex.free_symbols:
            return b**ex
        al = ex.coeff(L); be = sp.expand(ex - al*L)
        c = sp.nsimplify(b**al)
        out = b**be
        if c.is_negative:
            out *= SG; c = -c
        k = sp.log(c, 2)
        k = sp.nsimplify(k)
        assert k.is_rational, (b, ex, c)
        assert k.is_integer, (b, ex, c)
        return out*TW**int(k)
    def fix_exp(arg):
        arg = sp.expand(arg)
        if L not in arg.free_symbols:
            return sp.exp(arg)
        al = arg.coeff(L); rest = sp.expand(arg - al*L)
        q = sp.nsimplify(al/(sp.I*sp.pi))
        assert q.is_integer, arg
        return sp.exp(rest)*SG**int(q)
    e = e.replace(lambda x: x.is_Pow and L in x.exp.free_symbols, lambda x: fix_pow(x.base, x.exp))
    e = e.replace(lambda x: isinstance(x, sp.exp) and L in x.args[0].free_symbols, lambda x: fix_exp(x.args[0]))
    return e

def split_term(t):
    """t (a product) -> (key = transcendental monomial, rational coefficient in s, L, I)"""
    rat, tr = sp.Integer(1), sp.Integer(1)
    for f in sp.Mul.make_args(t):
        base, ex = f.as_base_exp()
        if f.is_number and not f.has(sp.pi):
            rat *= f
        elif f.free_symbols and f.free_symbols <= {S, L}:
            rat *= f
        else:
            if base == SG:
                f = SG**(int(ex) % 2)
            tr *= f
    return tr, rat

def exact_sum(terms, divisor):
    """sum_t t/divisor exactly: group by transcendental monomial (base Gammas, TW, SG, sqrt(pi), pi), sum the rational
    coefficients over a common denominator with polynomial arithmetic."""
    groups = {}
    dv = norm_powers(divisor)
    for t in terms:
        tt = sp.expand(sp.powsimp(norm_powers(t)/dv, force=True))
        for u in sp.Add.make_args(tt):
            k, r = split_term(u)
            groups.setdefault(k, []).append(r)
    out = 0
    for k, rats in groups.items():
        nds = [sp.fraction(sp.together(r)) for r in rats]
        den = sp.Integer(1)
        for _, d in nds:
            den = sp.lcm(den, sp.factor(d))
        num = sp.Integer(0)
        for n_, d in nds:
            num += sp.expand(n_*sp.cancel(den/d))
        out += k*sp.factor(sp.expand(num))/den
    return sp.simplify(out)

def ratio_exact(tot, divisor):
    """sum over the top-level terms of tot of term/divisor, reduced; returned as a cancelled rational function in
    (s, L, base Gammas, TW = 2^L, SG = (-1)^L, sqrt(pi)); SG^2 = 1 is imposed at the end."""
    RP = sp.Symbol('RP')
    e = sp.expand(tot, deep=False, power_exp=False, power_base=False, multinomial=False, log=False)
    dv = norm_powers(reduce_expr(divisor))
    acc = 0
    for t in sp.Add.make_args(e):
        r = sp.powsimp(norm_powers(reduce_expr(t))/dv, force=True)
        acc += r.subs(sp.sqrt(sp.pi), RP)
    out = sp.cancel(sp.together(acc))
    num, den = sp.fraction(out)
    num = sp.expand(num).subs(SG**2, 1); den = sp.expand(den).subs(SG**2, 1)
    return sp.factor(sp.cancel(num/den))

if __name__ == '__main__':
    import time
    res = {}
    for fn in ('extreme_pp.pkl', 'extreme_pm_odd_oddJ.pkl'):
        res.update(pickle.load(open(os.path.join(CACHE, fn), 'rb')))
    keys = sys.argv[1:] and [tuple(k.split(':')) for k in sys.argv[1:]] or sorted(res)
    for (name, which) in keys:
        t0 = time.time()
        Jsub = 2*L + 1 if name == 'oddJ' else 2*L
        tot = res[(name, which)].subs(J, Jsub)
        cl = claimed(name, which, Jsub)
        if cl == 0:
            terms = sp.Add.make_args(sp.expand(tot, deep=False, power_exp=False, power_base=False, multinomial=False, log=False))
            r = ratio_exact(tot, terms[0])
            print(f"{name:5s} {which}: absent sector: (amplitude)/(its first term) = {r}   [{time.time()-t0:.0f}s]", flush=True)
        else:
            r = ratio_exact(tot, cl)
            print(f"{name:5s} {which}: (symbolic-J amplitude)/(claimed closed form) = {r}   [{time.time()-t0:.0f}s]", flush=True)
