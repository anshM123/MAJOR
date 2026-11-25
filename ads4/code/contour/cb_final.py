"""
contour/cb_final.py (CONTOUR_BOUNDS.md §6): plug the certified contour numbers into §38's final inequality.

§38 structure (a >= 998, factors 100/a^2, 300/a^2 and |D(t + iY)| <= 5e12 kept as given):
  |E2| <= (100/a^2) M(L) + (4 sqrt w/Lam) X(L),
  X(L) = (300/a^2) [ e^{-3L/4} int_0^{1/4} |fB|(t + 3i/4) dt + int_{3/4}^{Y} |fB|(1/4 + iy) e^{-Ly} dy ] + 5e12 e^{-YL} int_{1/4}^oo |fB|(t + iY) dt,
  M(L) >= sqrt w (2/Lam) I_W(L),  I_W(L) = int_0^{3/4} W(y) e^{-Ly} dy      (A >= 0 dropped),
so  |E2|/M <= 100/a^2 + 2 X(L)/I_W(L) =: rho(L).
Certified inputs: cb_w.json (cell lower bounds of W), cb_short.json, cb_vert.json (cell sups of |f|, |B|, |fB|), cb_horiz_<Y>.json.
Variant 'sep'  : §38's factorisation, |fB| <= sup|f| sup|B| per piece (sup|B| over the piece), with the certified numbers.
Variant 'comb' : |fB| bounded jointly (int|fB| on the short and horizontal pieces, cell sups of |fB| on the vertical one).
Between grid points L_k < L < L_{k+1}: X is decreasing and I_W is decreasing in L, so rho(L) <= 100/a^2 + 2 X(L_k)/I_W(L_{k+1});
the thresholds below use this bound on a grid of step 0.01 (rigorous for L in the scanned range [5, 400]; beyond, every term of X
decays exponentially while I_W(L) >= W_lo(0)(1 - e^{-0.004 L})/L, so rho stays below its value at 400).
"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from flint import arb, ctx
ctx.prec = 128
CACHE = os.path.join(HERE, 'cache')
A = arb(998)
YS = [1.0, 1.25, 1.5, 1.75, 2.0, 2.25]


def ex_int(L, a, b):
    """int_a^b e^{-Ly} dy (exact, arb)."""
    return ((-L*a).exp() - (-L*b).exp())/L


class Data:
    def __init__(self):
        w = json.load(open(os.path.join(CACHE, 'cb_w.json')))
        c = w['cells']                                   # merge 10 cells (min of lower bounds) for speed
        self.wcells = [(arb(c[i][0]), arb(c[min(i + 9, len(c) - 1)][1]), arb(min(x[2] for x in c[i:i + 10])))
                       for i in range(0, len(c), 10)]
        self.w = w
        self.short = json.load(open(os.path.join(CACHE, 'cb_short.json')))['summary']
        v = json.load(open(os.path.join(CACHE, 'cb_vert.json')))
        self.vrows = v['rows']
        self.h = {}; self.hprof = {}
        for Y in YS:
            p = os.path.join(CACHE, f'cb_horiz_{Y}.json')
            if os.path.exists(p):
                d = json.load(open(p))
                self.h[Y] = d['summary']
                # D-profile variant: |D(t + iY)| <= 200 (1+t)^{5/4} for t <= 998 (below the turning point, g = 0), 5e12 beyond;
                # returned as the factor replacing 5e12 * int|fB|
                tot = arb(0)
                for r in d['rows']:
                    t1 = arb(r['t1'])
                    Dp = 200*(1 + t1)**(arb(5)/4) if r['t1'] <= 998 else arb('5e12')
                    tot += arb(r['intfB'])*Dp
                self.hprof[Y] = tot + arb('5e12')*arb(d['summary']['tail_intfB'])

    def IW(self, L):
        L = arb(L)
        return sum((lo*ex_int(L, a, b) for a, b, lo in self.wcells), arb(0))

    def pieces(self, L, Y, variant):
        L = arb(L)
        s = self.short
        vr = [r for r in self.vrows if r['y0'] < Y - 1e-12]
        if variant == 'sep':
            sB = max(s['supB'], max(r['supB'] for r in vr))   # one |B| bound on the detour, as §38
            short = arb(s['supf'])*arb(s['supB'])/4
            vert = sum((arb(r['supf'])*arb(sB)*ex_int(L, arb(r['y0']), arb(min(r['y1'], Y))) for r in vr), arb(0))
            hh = self.h[Y]
            hor = arb(hh['int_f'])*arb(hh['supB']) + arb(hh['tail_intf'])*arb(hh['tail_supB'])
        else:
            short = arb(s['int_fB'])
            vert = sum((arb(r['supfB'])*ex_int(L, arb(r['y0']), arb(min(r['y1'], Y))) for r in vr), arb(0))
            hor = arb(self.h[Y]['int_fB_total'])
        if variant == 'prof':
            P = dict(short=300/A**2*(-3*L/4).exp()*short, vert=300/A**2*vert, horiz=(-Y*L).exp()*self.hprof[Y])
        else:
            P = dict(short=300/A**2*(-3*L/4).exp()*short, vert=300/A**2*vert, horiz=arb('5e12')*(-Y*L).exp()*hor)
        return P

    def rho(self, Lk, Lk1, Y, variant):
        P = self.pieces(Lk, Y, variant)
        X = sum(P.values(), arb(0))
        return 100/A**2 + 2*X/self.IW(Lk1), P


def rho38(Lk, Lk1):
    """§38 as published (W >= 0.03 on [0, 0.01] only; |B| <= 20; crude sups; horizontal 4e31)."""
    Lk = arb(Lk); Lk1 = arb(Lk1)
    IW = arb('0.03')*(1 - (-Lk1/100).exp())/Lk1
    X = (-3*Lk/4).exp()/A**2*(arb('1.5e10') + arb('1.2e13')/(Lk - 17)) + arb('4e31')*(-9*Lk/4).exp()
    return 100/A**2 + 2*X/IW


def threshold(fn, target, Lmin=5, Lmax=400):
    """smallest L0 (to 0.01) such that the cell bound fn(L_k, L_{k+1}) < target on a cover of [L0, Lmax] by cells (step 0.5 from
    Lmax down until a cell fails, then step 0.01)."""
    L = Lmax
    while L - 0.5 >= Lmin and fn(L - 0.5, L) < target:
        L -= 0.5
    while L - 0.01 >= Lmin - 1e-9 and fn(L - 0.01, L) < target:
        L = round(L - 0.01, 2)
    return L


def main():
    D = Data()
    lines = []
    def pr(s=''):
        print(s, flush=True); lines.append(s)
    pr(__doc__)
    pr(f'I_W(L) lower bounds: ' + ', '.join(f'L={L}: {D.IW(L).lower()}' for L in (10, 20, 28, 40, 64)))
    pr(f'  (x L: ' + ', '.join(f'{float((D.IW(L)*L).lower()):.5f}' for L in (10, 20, 28, 40, 64)) + ')')
    res = {}
    for tgt in (0.1, 0.01, 0.001, 1.01e-4):
        f38 = lambda a, b: rho38(a, b)
        L38 = threshold(f38, tgt, Lmin=17.5)
        pr(f'§38 as published: rho < {tgt}: L >= {L38}')
        res[f'38_{tgt}'] = L38
    for Y in sorted(D.hprof):
        pr(f'horizontal leg Y = {Y}: 5e12*int|fB| <= {float((arb("5e12")*arb(D.h[Y]["int_fB_total"])).upper()):.4e};'
           f'  int|fB| x D-profile <= {float(D.hprof[Y].upper()):.4e}')
    for variant in ('sep', 'comb', 'prof'):
        for Y in sorted(D.h):
            if variant == 'prof' and Y != 2.25:
                continue
            for tgt in (0.1, 0.01, 0.001, 1.01e-4):
                L0 = threshold(lambda a, b: D.rho(a, b, Y, variant)[0], tgt, Lmin=2)
                r, P = D.rho(L0, L0 + 0.01, Y, variant)
                IW = D.IW(L0 + 0.01)
                share = {k: float((2*v/IW).upper()) for k, v in P.items()}
                pr(f'{variant:4s} Y = {Y:4}: rho < {tgt:8}: L >= {L0:7.2f};  at L0: rho <= {float(r.upper()):.3e}; '
                   f'contributions 2X/I_W: ' + ', '.join(f'{k} {v:.2e}' for k, v in share.items()) + f', 100/a^2 {float((100/A**2).upper()):.3e}')
                res[f'{variant}_{Y}_{tgt}'] = dict(L=L0, rho=float(r.upper()), share=share)
    for variant in ('sep', 'comb', 'prof'):
        for L in (20, 25, 28, 30, 40, 64):
            best = min((float(D.rho(L, L, Y, variant)[0].upper()), Y) for Y in D.h if variant != 'prof' or Y == 2.25)
            pr(f'{variant}: rho({L}) <= {best[0]:.4e} (best Y = {best[1]})')
            res[f'{variant}_rho_{L}'] = best
    json.dump(res, open(os.path.join(CACHE, 'cb_final.json'), 'w'), indent=1)
    open(os.path.join(CACHE, 'cb_final.txt'), 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
