"""
Audit of ANGULAR_ROUTE §21 item 3 (kernel bounds) and an extension of item 1 to J > 12.

(1) For J <= 12 the canonical structures/P_m come from code/kernels/data/structures_J*.json; for 14 <= J <= 24 they
    are recomputed with the repo pipeline (compute_kernels.py, cache redirected to ./cache) and normalised by
    canonical.canonical().
(2) Pipeline squared kernels |K|^2 = 12 |G P_m|^2 / N_{O,m} (m = 0 for pp/odd, m = 4 for pm/oddJ), with N from the
    independent pipeline.norm_poly (not the closed-form norm), compared with
      claim:  |K_even|^2 <~ 4^D 2^-J sqrt(J+1) sqrt(s h) / D^{15/2},   |K_odd|^2 <~ 4^D 2^-J sqrt(J+1)/(sqrt(s h) D^{15/2})
    at h = n + 1/4 (sin^2(pi h) = 1/2).
(3) Collinear identities for J = 14..24 (canonical structures from the pipeline).
(4) Closed-form scan of the ratio |K|^2 / claim over J <= 20000, 5 <= h <= 20000.
"""
import os, sys, json, pickle
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ['KERNELS_CACHE'] = os.path.join(HERE, 'cache')
sys.path.insert(0, os.path.join(HERE, '..', 'kernels'))
import sympy as sp, mpmath as mp
import kernels_lib as kl, closed_forms as cf, pipeline as pl
import ca1_collinear as c1
mp.mp.dps = 40
S = kl.S

def get_struct(J):
    """name -> (parity, monos, coeffs, K, {m: P_m})"""
    fn = os.path.join(HERE, '..', 'kernels', 'data', f'structures_J{J}.json')
    if os.path.exists(fn):
        d = json.load(open(fn)); out = {}
        for name, v in d.items():
            monos = [tuple(m[:6]) + (m[6],) for m in v['monomials']]
            out[name] = (v['parity'], monos, [sp.sympify(c, locals={'s': S}) for c in v['coefficients']],
                         tuple(v['K']), {int(m): sp.sympify(p, locals={'s': S}) for m, p in v['P'].items()})
        return out
    import canonical
    can = canonical.canonical(J)
    return {name: (par, monos, vec, K, P) for name, (par, monos, vec, K, P) in can.items()}

_norm_cache = {}
def norm_pipe(J, Delta, m):
    """pipeline N_{O,m} (norm_poly, symbolic in s), evaluated numerically"""
    if J not in _norm_cache:
        _norm_cache[J] = {mm: kl.poly_to_sym(v) for mm, v in pl.norm_poly(J).items()}
    nm = _norm_cache[J][m]
    sv = (mp.mpf(Delta) + J)/2
    val = mp.mpf(str(sp.N(sp.sympify(nm).subs(S, sp.Rational(str(sv))), 45)))
    return mp.pi**2.5*mp.mpf(2)**(4 - 2*Delta - 2*J)*val*mp.rgamma(Delta + J)*mp.rgamma(Delta + J - mp.mpf(1)/2)

def K2_pipe(J, struct, Delta, m):
    par, monos, vec, K, P = struct
    sv = (mp.mpf(Delta) + J)/2
    Pm = sp.sympify(P.get(m, 0))
    val = sp.N(Pm.subs(S, sp.Rational(str(sv))), 45)
    M = kl.G_prefactor(par, K, sv)*mp.mpc(str(sp.re(val)), str(sp.im(val)))
    return 12*abs(M)**2/norm_pipe(J, Delta, m)

def claim(J, Delta, odd):
    s = (Delta + J)/mp.mpf(2); h = (Delta - J)/mp.mpf(2)
    base = mp.mpf(4)**Delta*mp.mpf(2)**(-J)*mp.sqrt(J + 1)/Delta**mp.mpf(7.5)
    return base/mp.sqrt(s*h) if odd else base*mp.sqrt(s*h)

if __name__ == '__main__':
    Js = [int(x) for x in sys.argv[1:]] or [4, 6, 8, 10, 11, 12, 14, 15, 16, 18, 20, 24]
    print("(2)+(3): pipeline kernels vs closed forms vs claimed envelope; collinear identities")
    lpp, lpm, lod = sp.symbols('lpp lpm lodd')
    for J in Js:
        try:
            st = get_struct(J)
        except FileNotFoundError:
            print(f"J={J}: not computed"); continue
        # norm formula check
        dn = max(abs(norm_pipe(J, mp.mpf(J) + mp.mpf('7.3'), m)/cf.norm_O(mp.mpf(J) + mp.mpf('7.3'), J, m) - 1) for m in (0, 4, J))
        print(f"J={J:2d}: norm_poly vs closed-form norm, max rel dev {mp.nstr(dn, 3)}")
        for name, s_ in st.items():
            m0 = 4 if name in ('pm', 'oddJ') else 0
            cfun = {'pm': cf.Kpm2, 'pp': cf.Kpp2, 'odd': cf.Kodd2, 'oddJ': cf.KoddJ2}[name]
            odd = name in ('odd', 'oddJ')
            row = []
            for hh in (5.25, 10.25, 40.25, 200.25, 1000.25):
                D = mp.mpf(J) + 2*mp.mpf(hh)
                kp = K2_pipe(J, s_, D, m0)
                row.append((hh, kp/cfun(D, J), kp/claim(J, D, odd)))
            print(f"   {name:4s}: " + "  ".join(f"h={h}: pipe/cf={mp.nstr(a, 8)} pipe/claim={mp.nstr(b, 4)}" for h, a, b in row))
        if J >= 14:
            for opol in ('+',):
                Ft = {}
                for name, (par, monos, vec, K, P) in st.items():
                    lam = {'pp': lpp, 'pm': lpm, 'odd': lod, 'oddJ': lod}[name]
                    for k, v in c1.collinear_F(J, (par, monos, vec), opol).items():
                        Ft[k] = Ft.get(k, 0) + lam*v
                h = S - J; p = (h - 1)*(h - 2)*(h - 3)
                r22 = sp.factor(sp.simplify(Ft.get((-2, -2), 0)/(p*(h - 4))))
                r11 = sp.factor(sp.simplify(Ft.get((-1, -1), 0)/(-16*p)))
                r21 = sp.factor(sp.simplify(Ft.get((-2, -1), 0)/(-2*(2*h - 7)*(2*h - 5)*(2*h - 3))))
                print(f"   collinear: F22/[p(h)(h-4)] = {r22};  F11/[-16p] = {r11};  F21/[-2(2h-7)(2h-5)(2h-3)] = {r21}")
    print("(4) closed-form scan of |K|^2/claim  (h = n + 1/4)")
    for name, fun, odd, Jset in (('pp', cf.Kpp2, False, [4, 10, 100, 1000, 10000, 20000]),
                                 ('pm', cf.Kpm2, False, [6, 10, 100, 1000, 10000, 20000]),
                                 ('odd', cf.Kodd2, True, [4, 10, 100, 1000, 10000, 20000]),
                                 ('oddJ', cf.KoddJ2, True, [5, 11, 101, 1001, 10001, 20001])):
        lo, hi = mp.inf, 0
        for J in Jset:
            for hh in (5.25, 6.25, 10.25, 30.25, 100.25, 1000.25, 10000.25, 20000.25):
                D = mp.mpf(J) + 2*mp.mpf(hh)
                r = fun(D, J)/claim(J, D, odd)
                lo, hi = min(lo, r), max(hi, r)
        print(f"   {name:4s}: ratio in [{mp.nstr(lo, 4)}, {mp.nstr(hi, 4)}]")
