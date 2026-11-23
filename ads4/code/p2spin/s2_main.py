"""
p2spin/s2_main.py (P2_SPINNING.md §2): heavy symbols to O(1/m^2) of the conserved structures, helicity decomposition.

For q = e^{beta/2} in QS (exact rationals) and symbolic X, Y (rescaled x, y = X/m, Y/m), the polarisations are taken in the
helicity frame found in §2.1:  z1 = B(q^-1)(1, cos t1, sin t1),  z2' = B(q)(1, cos t2, sin t2)   (Cartesian (t, z, y), B(s) the boost
u -> s u, v -> v/s along n's spatial direction), and the symbol is expanded in E_i = e^{i t_i}:
    R_L(eps) = sum_{k = 0..2} eps^k sum_{(k1, k2)} R^(k)_{k1 k2}(X, Y) E1^k1 E2^k2 .
Outputs (cache/s2_main.pkl, cache/s2_main.txt):
  * leading amplitudes per mode (exact in q, from a symbolic-q run at numeric X, Y);
  * helicity basis (even J): 'same' (modes (2,2),(-2,-2)) and 'opp' (modes (2,-2),(-2,2)), unit leading amplitude;
  * for each helicity structure and order k = 1, 2: the polynomial R^(k)_{k1 k2}(X, Y)/R^(0)_{lead} per mode, split into
    (i) leading-mode constant (congruence/normalisation), (ii) leading-mode X,Y-dependent, (iii) other modes (helicity leakage);
  * the reflection identity R_R(eps) = R_L(-eps) checked numerically on the actual right-side construction.
"""
import os, sys, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sympy as sp
import p2_series as PS
PS.ORDER = 2
import s2_symbols as S2

E1, E2 = sp.symbols('E1 E2')
QS = [sp.Rational(11, 10), sp.Rational(3, 2), sp.Rational(2)]
OUT = []
pr = lambda *a: (print(' '.join(str(x) for x in a), flush=True), OUT.append(' '.join(str(x) for x in a)))


def hel_pol(E, s, qq):
    return (qq**s*(1 + (E + 1/E)/2), qq**(-s)*(1 - (E + 1/E)/2), (E - 1/E)/(2*sp.I))


A_ = sp.symbols('pa0:3'); B_ = sp.symbols('pb0:3')


def _mode_table(qv):
    """map from the pol monomials a_i a_j b_k b_l to helicity-frame modes (k1, k2)."""
    h1 = hel_pol(E1, -1, qv); h2 = hel_pol(E2, 1, qv)
    tab = {}
    for i in range(3):
        for j in range(i, 3):
            for k in range(3):
                for l in range(k, 3):
                    g = sp.expand(h1[i]*h1[j]*h2[k]*h2[l]*E1**2*E2**2)
                    P = sp.Poly(g, E1, E2)
                    tab[(i, j, k, l)] = {(a - 2, b - 2): c for (a, b), c in P.terms()}
    return tab


_TABS = {}


def modes_from_pol(expr, qv):
    """expr: polynomial in A_, B_ (degree (2, 2)) with coefficients in X, Y -> {mode: coefficient}."""
    if qv not in _TABS:
        _TABS[qv] = _mode_table(qv)
    tab = _TABS[qv]
    P = sp.Poly(sp.expand(expr), *A_, *B_)
    out = {}
    for mon, c in P.terms():
        ia = [i for i in range(3) for _ in range(mon[i])]; ib = [k for k in range(3) for _ in range(mon[3 + k])]
        assert len(ia) == 2 and len(ib) == 2
        for md, w in tab[(ia[0], ia[1], ib[0], ib[1])].items():
            out[md] = out.get(md, 0) + w*c
    return {k: sp.expand(v) for k, v in out.items() if sp.expand(v) != 0}


def sym_modes(side, key, cf, xv, yv, qv):
    s = S2.symbol(side, key[0], cf, xv, yv, A_, B_, qv)
    return s.v, [modes_from_pol(c, qv) for c in s.c]


def main():
    t0 = time.time()
    can = S2.CAN()
    fams = [('even', 1, 0), ('even', 1, 1), ('odd', 1, 0), ('odd', -1, 0)]
    res = {}
    # ---- leading amplitudes, exact in q (numeric X, Y; leading is X, Y independent -- checked at two points)
    qs = sp.Symbol('q', positive=True)
    P1 = ((sp.Rational(3, 2), sp.Rational(2, 3), sp.Rational(1, 3)), (1, sp.Rational(5, 4), sp.Rational(-1, 2)))
    P2 = ((sp.Rational(1, 2), sp.Rational(7, 3), sp.Rational(-2, 3)), (3, sp.Rational(1, 4), sp.Rational(1, 5)))
    PS.ORDER = 0
    lead = {}
    for sec, sig, k in fams:
        cf = can[(sec, sig)][k]
        v1, m1 = sym_modes('L', (sec, sig), cf, *P1, qs); v2, m2 = sym_modes('L', (sec, sig), cf, *P2, qs)
        same = all(sp.simplify(m1[0].get(kk, 0) - m2[0].get(kk, 0)) == 0 for kk in set(m1[0]) | set(m2[0]))
        lead[(sec, sig, k)] = {kk: sp.factor(c) for kk, c in m1[0].items()}
        pr(f'[{sec} sigma={sig:+d} basis {k}] valuation {v1}; leading X,Y-independent: {same}; modes: '
           + ', '.join(f'{kk}: {sp.factor(c)}' for kk, c in sorted(m1[0].items())))
    # helicity basis for even J
    a_s = lead[('even', 1, 0)][(2, 2)]; a_o = lead[('even', 1, 0)][(2, -2)]
    b_s = lead[('even', 1, 1)][(2, 2)]; b_o = lead[('even', 1, 1)][(2, -2)]
    det = sp.factor(a_s*b_o - a_o*b_s)
    # same = (b_o*basis0 - a_o*basis1)/det ; opp = (-b_s*basis0 + a_s*basis1)/det
    HB = {'same': (sp.factor(b_o/det), sp.factor(-a_o/det)), 'opp': (sp.factor(-b_s/det), sp.factor(a_s/det))}
    pr(f'helicity basis (even J): same = {HB["same"][0]}*b0 + {HB["same"][1]}*b1 ;  opp = {HB["opp"][0]}*b0 + {HB["opp"][1]}*b1')
    res['lead'] = lead; res['HB'] = HB
    pr(f'   [{time.time() - t0:.0f}s]')
    # ---- orders 1, 2, symbolic X, Y, at q in QS
    PS.ORDER = 2
    X = S2.xs; Y = S2.ys
    res['orders'] = {}
    for qv in QS:
        sm = {}
        for sec, sig, k in fams:
            v, ms = sym_modes('L', (sec, sig), can[(sec, sig)][k], X, Y, qv)
            sm[(sec, sig, k)] = (v, ms)
        # helicity combos
        hb = {}
        for h in ('same', 'opp'):
            c0, c1 = [c.subs(qs, qv) for c in HB[h]]
            comb = []
            for o in range(3):
                d = {}
                for kk in set(sm[('even', 1, 0)][1][o]) | set(sm[('even', 1, 1)][1][o]):
                    d[kk] = sp.expand(c0*sm[('even', 1, 0)][1][o].get(kk, 0) + c1*sm[('even', 1, 1)][1][o].get(kk, 0))
                comb.append({kk: c for kk, c in d.items() if c != 0})
            hb[h] = comb
        hb['OO'] = sm[('odd', 1, 0)][1]; hb['odd'] = sm[('odd', -1, 0)][1]
        res['orders'][qv] = hb
        pr(f'q = {qv} (beta = {float(2*sp.log(qv)):.4f})   [{time.time() - t0:.0f}s]')
        for h, comb in hb.items():
            L0 = comb[0]
            lead_modes = sorted(L0)
            ref = L0[lead_modes[-1]]
            pr(f'   {h:5s}: leading modes {lead_modes} (amp {[sp.nsimplify(L0[kk]/ref) for kk in lead_modes]})')
            for o in (1, 2):
                Lo = comb[o]
                cst = {kk: Lo[kk].subs({**{x: 0 for x in X}, **{y: 0 for y in Y}}) for kk in Lo}
                lm = {kk: sp.expand(Lo[kk]/L0[kk]) for kk in Lo if kk in L0}
                other = {kk: Lo[kk] for kk in Lo if kk not in L0}
                pr(f'      order {o}: leading-mode relative corrections:')
                for kk, c in lm.items():
                    c0 = c.subs({**{x: 0 for x in X}, **{y: 0 for y in Y}})
                    pr(f'         {kk}: const {sp.nsimplify(c0)}; X,Y part: {sp.expand(c - c0)}')
                pr(f'      order {o}: leakage modes {sorted(other)}; max |const| / |lead| = '
                   f'{max([abs(complex(sp.N(cst[kk]/ref))) for kk in other] or [0]):.3g}')
    # ---- reflection identity R_R(eps) = R_L(-eps) (numeric check of the construction)
    PS.ORDER = 2
    for sec, sig, k in fams:
        cf = can[(sec, sig)][k]
        sl = S2.symbol('L', sec, cf, *P1, hel_pol(E1, -1, QS[1]), hel_pol(E2, 1, QS[1]), QS[1])
        sr = S2.symbol('R', sec, cf, *P1, hel_pol(E1, -1, QS[1]), hel_pol(E2, 1, QS[1]), QS[1])
        ok = sl.v == sr.v and all(sp.simplify(sr.c[o] - (-1)**(o + sl.v - sr.v)*sl.c[o]*(-1)**sl.v) == 0 for o in range(3))
        pr(f'reflection R_R(eps) = R_L(-eps) [{sec} sigma={sig:+d} basis {k}]: {ok}')
    pickle.dump({'lead': {k: {kk: sp.srepr(c) for kk, c in v.items()} for k, v in res['lead'].items()},
                 'HB': {h: tuple(sp.srepr(c) for c in v) for h, v in res['HB'].items()},
                 'orders': {sp.srepr(qv): {h: [{kk: sp.srepr(c) for kk, c in d.items()} for d in comb] for h, comb in hb.items()}
                            for qv, hb in res['orders'].items()}},
                open(os.path.join(HERE, 'cache', 's2_main.pkl'), 'wb'))
    open(os.path.join(HERE, 'cache', 's2_main.txt'), 'w').write('\n'.join(OUT) + '\n')


if __name__ == '__main__':
    main()
