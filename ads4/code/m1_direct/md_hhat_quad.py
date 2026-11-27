"""md_hhat_quad.py -- independent evaluation of the k = 2 stripped light action  hhat^(2)_b = T * int g(nu) S_b(nu) dnu  by direct
nu-quadrature of the PUBLISHED continuum weight (prl/data/k2/k2_continuum_D_<Dg>_rep.json, via m1_remainder/m1_weight.py, read-only)
against the ladder polynomial S_b (md_ladder, own sympy implementation) -- i.e. without the certified moment balls of fullpos.
g(nu) = -(9216/Lam) (4nu^2+1)/((4nu^2+9)(4nu^2+49)) chi(nu) theta(nu)   [Eq. (fp-moments)].
Compares with m1_remainder/no_eft per-state values (moments route, k = 2 + k = 3).
Usage: python3 md_hhat_quad.py Dg  -> out/hhat_quad_<Dg>.txt
"""
import os, sys, math, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'm1_remainder'))
import numpy as np
import mpmath as mp
import sympy as sp
import md_ladder as ML
import m1_weight as MW

mp.mp.dps = 60


def main(Dg):
    W = MW.Weight(Dg)
    Lam = mp.sqrt(W.lam)
    numax = float(Lam) + 170
    # Gauss-Legendre panels of width 2 (the weight's Gaussian basis has widths >= ~0.45 near the top; checked by panel doubling)
    xg, wg = np.polynomial.legendre.leggauss(24)
    edges = np.arange(0.0, numax + 2, 2.0)
    nodes, wts = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        nodes += list(0.5*(b - a)*xg + 0.5*(a + b)); wts += list(0.5*(b - a)*wg)
    t0 = time.time()
    gv = []
    for v in nodes:
        nu = mp.mpf(v)
        gv.append(-(mp.mpf(9216)/Lam)*(4*nu**2 + 1)/((4*nu**2 + 9)*(4*nu**2 + 49))*W.chi_mp(nu)*W.theta_mp(nu))
    print(f'weight at {len(nodes)} nodes: {time.time() - t0:.0f}s', flush=True)
    Z, Y = ML.zy(60)
    L = [f'# hhat^(2)_b = T int g S_b dnu by direct quadrature of the published Delta_g = {Dg} weight; GL24 on width-2 panels to Lam+170']
    for (tw, J, b) in [(2, 8, 0), (2, 8, 4), (2, 20, 0), (2, 20, 4), (3, 8, 0), (4, 8, 0), (4, 20, 0), (2, 4, 0), (998, 8, 0), (1000, 8, 0)]:
        if J - b < 0:
            continue
        Delta = tw + J; T = (Delta - 1)*(Delta - 2)
        S = sp.Poly(ML.h2(b, Delta, J, Z, Y), ML.s)
        co = [mp.mpf(sp.Rational(c).p)/mp.mpf(sp.Rational(c).q) for c in S.all_coeffs()]
        tot = mp.mpf(0); ab = mp.mpf(0)
        for v, w, g in zip(nodes, wts, gv):
            s = mp.mpf(v)**2 + mp.mpf(1)/4
            val = mp.polyval(co, s)
            tot += w*g*val; ab += abs(w*g*val)
        hh = T*tot
        L.append(f'  twist {tw:4d} J {J:3d} b {b}: hhat^(2) = {mp.nstr(hh, 12)}  log10|.| = {float(mp.log10(abs(hh))):.3f}   '
                 f'cancellation int|gS|/|int gS| = {mp.nstr(ab/abs(tot), 3)}')
        print(L[-1], flush=True)
    open(os.path.join(HERE, 'out', f'hhat_quad_{Dg}.txt'), 'w').write('\n'.join(L) + '\n')


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 1000)
