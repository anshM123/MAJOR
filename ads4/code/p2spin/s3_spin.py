"""
p2spin/s3_spin.py (P2_SPINNING.md §3.2): O(1/m^2) ring-class coefficients of the spinning heavy image in the
HELICITY-DIAGONAL KERNEL MODEL [M], via the CHMRS mechanism validated in s3_mech (exact).

Model [M]: the EE kernel acts on each helicity-mode class t in {same, opp} (modes (+-2, +-2) of s2_main's helicity frame) as the
scalar Regge-moment kernel times a mode-class normalisation, and annihilates every other mode ('leakage' modes, reported separately).
At leading order this reproduces the exact facts used by the project (pure-helicity leading symbols, s2_main; helicity-diagonal
leading image, EE_FUNCTIONAL §0 item 3).  Beyond leading order it is NOT derived: the kernel has three polarisation channels
(LARGESPIN_W §9.1) and the model keeps one.

For structures h, h' and mode class t, with normalised left coefficients L_{h,t}(eps) = R_{h,t}/R^(0)_{t,t} (s2_main; average of the
two modes of t, which are complex conjugates) and the reflection identity R_R(eps) = R_L(-eps) (s2_main, checked):
    M_{hh'}  ~  sum_t n_t [L_{h,t}(eps) L_{h',t}(-eps)]  x  seed (1 + eps^2 s2)
The eps^2 part is a polynomial in X, Y; the mechanism maps it to c0 P + c1 P' at the state's eta = cosh beta.
Output: per q, the nu^2 coefficients of P and P' (classes B and C2) and the nu^0 P' coefficient, for the seed (= scalar CHMRS,
validated) and for the spinning symbol, diagonal and off-diagonal; Lambda_amp(model) as operator norms over the 2x2 (same, opp)
block; the leakage (other modes) sizes.  -> cache/s3_spin.txt, cache/s3_spin.pkl
"""
import os, sys, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sympy as sp
import numpy as np
import s3_mech as M

OUT = []
pr = lambda *a: (print(' '.join(str(x) for x in a), flush=True), OUT.append(' '.join(str(x) for x in a)))
xs = sp.symbols('xu xv xy'); ys = sp.symbols('yu yv yy')
LC = {xs[0]: M.X[0] + M.X[1], xs[1]: M.X[0] - M.X[1], xs[2]: M.X[2],
      ys[0]: M.Y[0] + M.Y[1], ys[1]: M.Y[0] - M.Y[1], ys[2]: M.Y[2]}
CLASSES = {'same': [(2, 2), (-2, -2)], 'opp': [(2, -2), (-2, 2)]}
A = 7       # Fourier-Laplace power for k = 2, Delta_phi = 3 (validated: alpha = 2 Delta_phi + k - 1)


def load():
    d = pickle.load(open(os.path.join(HERE, 'cache', 's2_main.pkl'), 'rb'))
    ev = lambda s: sp.sympify(eval(s, {**sp.__dict__}))
    orders = {ev(qk): {h: [{tuple(kk): ev(c) for kk, c in dd.items()} for dd in comb] for h, comb in hb.items()}
              for qk, hb in d['orders'].items()}
    return orders


def ring(poly, P, Pb):
    """(nu^2 coeff of P, nu^0 coeff of P, nu^2 coeff of P', nu^0 coeff of P') of the mechanism image of poly."""
    if poly == 0:
        return (0, 0, 0, 0)
    expr, et0, reps = M.apply_insertion(sp.expand(poly), P, Pb, A)
    c0, c1 = M.reduce_legendre(expr, reps, et0)
    c0 = sp.expand(c0); c1 = sp.expand(c1)
    return (c0.coeff(M.nu2, 1), c0.coeff(M.nu2, 0), c1.coeff(M.nu2, 1), c1.coeff(M.nu2, 0))


def series_cls(comb, t):
    """normalised class-t coefficient series [L0, L1, L2] (average over the two complex-conjugate modes)."""
    lead_amp = None
    out = []
    for o in range(3):
        vals = [comb[o].get(md, 0) for md in CLASSES[t]]
        out.append(vals)
    return out


def main():
    t0 = time.time()
    orders = load()
    res = {}
    for qv, hb in orders.items():
        P, Pb, (s1, s2) = M.seed_series(qv)
        eta0 = sp.nsimplify((qv**2 + qv**-2)/2)
        pr(f'q = {qv}: beta = {float(2*sp.log(qv)):.4f}, eta = {eta0}   [{time.time() - t0:.0f}s]')
        seed = ring(s2, P, Pb)
        pr(f'   seed (scalar CHMRS, Regge moment): nu^2 P: {sp.nsimplify(seed[0])} (= {float(seed[0]/eta0):.6f} eta), '
           f'nu^2 P\': {sp.nsimplify(seed[2])} (= {float(seed[2]/(eta0**2 - 1)):.6f} (eta^2-1)), nu^0 P\': {float(seed[3]):.6f}')
        # leading amplitudes of each class from the pure-helicity structures
        ref = {t: hb[t][0][CLASSES[t][0]] for t in ('same', 'opp')}
        # left series per (structure h, class t): average over the two conjugate modes, normalised to class-t leading amplitude
        Ls = {}
        for h in ('same', 'opp'):
            for t in ('same', 'opp'):
                ser = []
                for o in range(3):
                    v = [hb[h][o].get(md, 0)/hb[t][0][md] for md in CLASSES[t]]
                    ser.append([sp.expand(sp.expand(x).subs(LC)) for x in v])
                Ls[(h, t)] = ser
        mat = {}
        for pairing in ('direct', 'conj'):
          for h in ('same', 'opp'):
            for h2 in ('same', 'opp'):
                for t in ('same', 'opp'):
                    tot = 0; o1 = 0
                    for im in range(2):                       # each of the two conjugate modes, averaged
                        jm = im if pairing == 'direct' else 1 - im     # right mode paired with left mode im
                        l = [Ls[(h, t)][o][im] for o in range(3)]; r = [Ls[(h2, t)][o][jm]*(-1)**o for o in range(3)]
                        tot += (l[0]*r[2] + l[1]*r[1] + l[2]*r[0])/2
                        o1 += (l[0]*r[1] + l[1]*r[0])/2
                    tot = sp.expand(tot); o1 = sp.expand(o1)
                    re = sp.expand((tot + tot.subs(sp.I, -sp.I))/2)          # X, Y real
                    im_part = sp.expand(tot - re)
                    zero = {v: 0 for v in M.X + M.Y}
                    cst = re.subs(zero)                      # X,Y-independent part: basis/normalisation drift (congruence class, W1c)
                    rr = list(ring(sp.expand(re - cst), P, Pb)); rr.append(cst)
                    mat[(pairing, h, h2, t)] = dict(ring=rr, o1=o1, im=im_part, lead=sp.expand((lambda z: (z + z.subs(sp.I, -sp.I))/2)(sp.expand(Ls[(h, t)][0][0]*Ls[(h2, t)][0][0]))))
                    pr(f'   [{pairing}] M[{h},{h2}] via class {t}: leading weight {sp.nsimplify(mat[(pairing, h, h2, t)]["lead"])}; '
                       f'O(eps) part {"0" if o1 == 0 else str(sp.factor(o1))[:70]}; eps^2 symbol -> nu^2 P: {float(rr[0]):+.5f}, '
                       f'nu^0 P: {float(rr[1]):+.5f}, nu^2 P\': {float(rr[2]):+.5f}, nu^0 P\': {float(rr[3]):+.5f}'
                       + ('' if im_part == 0 else f'  [imag part: {str(im_part)[:50]}]'))
        # leakage: modes outside the |h| = 2 classes, at orders 1, 2 (constant term size and X,Y polynomial size)
        for h in ('same', 'opp'):
            for o in (1, 2):
                lk = {md: c for md, c in hb[h][o].items() if md not in CLASSES['same'] + CLASSES['opp']}
                mx = 0
                for md, c in lk.items():
                    cc = sp.Poly(sp.expand(c/ref[h]), *xs, *ys)
                    mx = max(mx, max(abs(complex(sp.N(v))) for v in cc.coeffs()))
                pr(f'   leakage of {h} at order {o}: modes {sorted(lk)}; max |polynomial coefficient| / leading = {mx:.3g}')
        res[qv] = dict(seed=seed, mat=mat, eta=eta0)
        # Lambda_amp(model): spinning total = seed (x leading weight) + symbol, per class; operator norm over (same, opp)
        for pairing in ('direct', 'conj'):
            for k, lab in ((0, 'B  (nu^2 P)'), (2, 'C2 (nu^2 P\')'), (1, 'nu^0 P [X,Y-dep.]'), (3, 'nu^0 P\''), (4, 'const (normalisation)')):
                Mx = np.zeros((2, 2)); Ms = np.zeros((2, 2))
                for i, h in enumerate(('same', 'opp')):
                    for j, h2 in enumerate(('same', 'opp')):
                        Ms[i, j] = float(sum(mat[(pairing, h, h2, t)]['ring'][k] for t in ('same', 'opp')))
                        Mx[i, j] = Ms[i, j] + (float(seed[k]) if k < 4 else 0.0)*float(sum(mat[(pairing, h, h2, t)]['lead'] for t in ('same', 'opp')))
                ev = np.linalg.eigvalsh((Mx + Mx.T)/2); evs = np.linalg.eigvalsh((Ms + Ms.T)/2)
                pr(f'   [{pairing}] {lab}: symbol-only 2x2 = {np.round(Ms, 5).tolist()} ; seed+symbol eigenvalues {np.round(ev, 5).tolist()}; '
                   f'scalar {float(seed[k]) if k < 4 else float("nan"):.5f}; ratio max|eig|/|scalar| = {max(abs(ev))/abs(float(seed[k])) if k < 4 else float("nan"):.3f}; '
                   f'symbol-only max|eig|/|scalar| = {max(abs(evs))/abs(float(seed[k])) if k < 4 else float("nan"):.3f}')
    pickle.dump({sp.srepr(q): {'seed': [sp.srepr(sp.sympify(x)) for x in r['seed']],
                               'mat': {k: [sp.srepr(sp.sympify(x)) for x in v['ring']] for k, v in r['mat'].items()},
                               'lead': {k: sp.srepr(sp.sympify(v['lead'])) for k, v in r['mat'].items()}}
                 for q, r in res.items()}, open(os.path.join(HERE, 'cache', 's3_spin.pkl'), 'wb'))
    open(os.path.join(HERE, 'cache', 's3_spin.txt'), 'w').write('\n'.join(OUT) + '\n')


if __name__ == '__main__':
    main()
