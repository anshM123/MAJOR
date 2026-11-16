"""
Parse the CHMRS ancillary file regge_moments_order8.m (arXiv:2106.10274, anc/) into sympy.

Defines (CHMRS comments, verbatim meaning):
  seriesPi : Regge expansion of a_{Df} m^{2k}/gamma_{2Df+k-1}(nu)^2 * Pihat_{k,nu}[Delta,J]  (heavy density, i.e. with the
             positive factor q_{Delta,J} and 2 sin^2 removed), to order 1/m2^4
  seriesC  : Regge expansion of m^{2k} C_{k,nu}[Delta,J], to order 1/m2^4
  seriesCfixedJ : the same at fixed J (bulk-point regime), polynomial in nu^2 and J2 = J(J+d-2)
Variables: m2 = (Delta-J-d+1)(Delta+J-1), etaAdS = cosh beta, J2 = J(J+d-2); pj[j,x,n] = n-th x-derivative of the
Gegenbauer harmonic function normalised to 1 at x = 1 (d = 3: Legendre P_{-1/2+i nu}).
The ancillary file is copied into cache/ (ignored); this module only reads it.
"""
import re, os, pickle, sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, 'cache', 'regge_moments_order8.m')
PKL = os.path.join(HERE, 'cache', 'chmrs_series.pkl')

d, k, Df, al, nu, eta, m2, J2 = sp.symbols('d k Df al nu eta m2 J2')
P = sp.symbols('P0:9')          # P[n] = (d/deta)^n P_{(2-d)/2+i nu}(eta)


def _convert(txt):
    txt = txt.replace('\\[CapitalDelta]\\[Phi]', 'Df').replace('\\[Alpha]', 'al').replace('\\[Nu]', 'nu')
    txt = txt.replace('\\[Eta]AdS', 'eta')
    txt = re.sub(r'pj\[\(2 - d\)/2 \+ I\*nu,\s*eta,\s*(\d+)\]', r'P\1', txt)
    txt = txt.replace('^', '**').replace('\n', ' ')
    return txt


def load():
    if os.path.exists(PKL):
        return pickle.load(open(PKL, 'rb'))
    raw = open(SRC).read()
    # split the three definitions
    keys = (('seriesPi', 'series\\[CapitalPi] ='), ('seriesC', 'seriesC ='), ('seriesCfixedJ', 'seriesCfixedJ ='))
    pos = [raw.index(key) for _, key in keys] + [len(raw)]
    parts = {name: raw[pos[i] + len(key):pos[i + 1]] for i, (name, key) in enumerate(keys)}
    loc = {'d': d, 'k': k, 'Df': Df, 'al': al, 'nu': nu, 'eta': eta, 'm2': m2, 'J2': J2, 'I': sp.I}
    for n in range(9):
        loc['P%d' % n] = P[n]
    out = {}
    for name, txt in parts.items():
        expr = sp.sympify(_convert(txt), locals=loc)
        # organise by powers of 1/m2
        expr = sp.expand(expr.subs(m2, 1/sp.Symbol('eps')))
        eps = sp.Symbol('eps')
        out[name] = [sp.expand(expr.coeff(eps, j)) for j in range(5)]
    pickle.dump(out, open(PKL, 'wb'))
    return out


if __name__ == '__main__':
    import time
    t0 = time.time()
    s = load()
    for name, cs in s.items():
        print(name, [len(sp.Add.make_args(c)) for c in cs], 'free:', sorted(str(x) for x in set().union(*[c.free_symbols for c in cs])))
    print('time', time.time() - t0)
