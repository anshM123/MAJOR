"""p4win/p4w_series.py (P4_WINDOW.md §2): CHMRS's scalar Regge-moment series ('series Pi' of the ancillary regge_moments_order8.m,
Delta_phi = 3, k = 2, alpha = 7, d = 3) to ALL its orders 1/m^2 ... 1/m^8, per unit leading image (divided by 2, so the leading term is
P = P_{-1/2+i nu}(eta)), Legendre-reduced to  T_k = r0_k(eta, nu^2) P + r1_k(eta, nu^2) P'(eta)  (P' = dP/deta), and the flat-scaling
part g_k(b) of each order (P -> J0(b), P' -> -(nu^2/b) J1(b), eta^2 - 1 -> b^2/nu^2; the terms with nu^(2k) net scaling).
Read-only: the ancillary (scratch copy, argv[1]); p2con/pc_flat.reduce and flat_part logic re-implemented here (orders > 2 need more
terms of the w-series than pc_flat keeps).
-> cache/series.pkl  (sympy srepr strings), cache/series.txt
Usage: python3 p4w_series.py <regge_moments_order8.m>
"""
import os, sys, re, pickle, time
sys.dont_write_bytecode = True
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
nu2, eta, w, bsym, J0s, J1s = sp.symbols('nu2 eta w b J0 J1')
OUT = []
pr = lambda s: (print(s, flush=True), OUT.append(s))


def parse(path, d=3, Dphi=3, alpha=7):
    txt = open(path).read()
    s = txt[txt.index('series\\[CapitalPi] ='):txt.index('seriesC =')]
    s = s[s.index('=') + 1:].strip().rstrip(';')
    for a, b in (('\\[Eta]AdS', 'eta'), ('\\[Nu]', 'nu'), ('\\[CapitalDelta]\\[Phi]', 'Dphi'), ('\\[Alpha]', 'alpha'), ('^', '**')):
        s = s.replace(a, b)
    s = re.sub(r'\s+', ' ', s)
    s = re.sub(r'pj\[\s*\(2 - d\)/2 \+ I\*nu,\s*eta,\s*(\d+)\s*\]', r'G\1', s)
    assert 'pj' not in s
    nu, m2 = sp.symbols('nu m2')
    loc = {'nu': nu, 'eta': eta, 'm2': m2, 'd': sp.Integer(d), 'Dphi': sp.Integer(Dphi), 'alpha': sp.Integer(alpha), 'I': sp.I}
    loc.update({f'G{k}': sp.Symbol(f'G{k}') for k in range(13)})
    ex = sp.expand(sp.sympify(s, locals=loc))
    free = ex.free_symbols - {nu, eta, m2} - {sp.Symbol(f'G{k}') for k in range(13)}
    assert not free, free
    terms = {}
    for k in range(0, 5):
        t = sp.expand(ex.coeff(m2, -k)/2)
        terms[k] = sp.expand(t.subs(nu, sp.sqrt(nu2)))
    return terms


def reduce(t):
    """G_n = d^n P/d eta^n -> (r0 P + r1 P') via Legendre's equation (eta^2-1)P'' + 2 eta P' + (nu^2 + 1/4) P = 0."""
    sol = {0: (sp.Integer(1), sp.Integer(0)), 1: (sp.Integer(0), sp.Integer(1))}
    for k in range(0, 11):
        c0 = -(2*(k + 1)*eta*sol[k + 1][0] + (k*(k + 1) + nu2 + sp.Rational(1, 4))*sol[k][0])/(eta**2 - 1)
        c1 = -(2*(k + 1)*eta*sol[k + 1][1] + (k*(k + 1) + nu2 + sp.Rational(1, 4))*sol[k][1])/(eta**2 - 1)
        sol[k + 2] = (sp.cancel(c0), sp.cancel(c1))
    r0 = sp.cancel(sum(t.coeff(sp.Symbol(f'G{k}'))*sol[k][0] for k in range(13)))
    r1 = sp.cancel(sum(t.coeff(sp.Symbol(f'G{k}'))*sol[k][1] for k in range(13)))
    return r0, r1


def flat_part(r0, r1, order, nw=14):
    """flat-scaling terms at O(y^order): nu^{2j} w^i P with j - i = order; nu^{2j} w^i P' with j - i + 1 = order (w = eta^2 - 1)."""
    tot = 0
    for r, kind in ((r0, 'P'), (r1, 'D')):
        e = sp.expand(sp.series(sp.cancel(r.subs(eta, sp.sqrt(1 + w))), w, 0, nw).removeO())
        for term in sp.Add.make_args(e):
            pw = term.as_powers_dict(); i = int(pw.get(w, 0)); jn = int(pw.get(nu2, 0))
            c = term/(w**i*nu2**jn)
            assert not c.free_symbols, term
            if kind == 'P' and jn - i == order:
                tot += c*bsym**(2*i)*J0s
            if kind == 'D' and jn - i + 1 == order:
                tot += -c*bsym**(2*i - 1)*J1s
            if kind == 'P' and jn - i > order or kind == 'D' and jn - i + 1 > order:
                raise ValueError('super-flat term', term, kind)
    return sp.expand(tot)


def main():
    t0 = time.time()
    terms = parse(sys.argv[1])
    assert sp.simplify(terms[0] - sp.Symbol('G0')) == 0
    res = {}
    for k in range(1, 5):
        r0, r1 = reduce(terms[k])
        g = flat_part(r0, r1, k)
        res[k] = dict(r0=sp.srepr(r0), r1=sp.srepr(r1), g=sp.srepr(g))
        pr(f'order 1/m^{2*k}: nu-degree of r0 {sp.Poly(sp.numer(sp.together(r0)), nu2).degree()}, flat part g_{k}(b) = {g}   [{time.time() - t0:.0f}s]')
    pickle.dump(res, open(os.path.join(HERE, 'cache', 'series.pkl'), 'wb'))
    open(os.path.join(HERE, 'cache', 'series.txt'), 'w').write('\n'.join(OUT) + '\n')


if __name__ == '__main__':
    main()
