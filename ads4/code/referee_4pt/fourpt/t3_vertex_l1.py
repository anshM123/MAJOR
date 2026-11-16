"""
Heavy limit of the doubly-spinning shadow integrand for <JJJJ> (l = 1), pointwise at fixed (X, Y) = (m x, m y).
For the Weyl-covariant basis (v_1, v_2) of t1 (right structure = (1/b^(0)) x the same family at (3-Delta,-1-J)):
   T^{ab}(z1,z2',z3,z4') / [((-1)^J/(4 b^(0))) 2 4^4 |x|^4 |y|^4 e^{2P.x + 2Pbar.y}]  ->  sum_h V_a(h1,h2) V_b(h3,h4) prod (z.eps_h)
with the closed forms of predict4.V_l1.  Checks: (i) ratio exact/closed form for all 16 helicity components and all (a,b) -> 1 as
1 + O(1/m); (ii) longitudinal leakage -> 0; (iii) the O(1/m) term is antisymmetric under the exchange of the two vertices'
helicity labels together with (a <-> b) (Weyl reflection), so the symmetrised combination approaches 1 like 1/m^2.
Also: (iv) checks of the kernel conventions: K_{mu,0} = F^nu_mu (e2e one-pair kernel).
"""
import sys, pickle, os
import mpmath as mp, sympy as sp
from vertex2 import helicity_components
from predict4 import V_l1, weyl
from moments4 import kernel_mp, kernel_np
from embed2 import Delta, J, CACHE
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'e2e'))
from formulas import Fnu_mp                   # noqa: E402
mp.mp.dps = 40
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)

ws, vecs = pickle.load(open(os.path.join(CACHE, 'consbasis_l1.pkl'), 'rb'))
def coeffs(v, D, Jv):
    return {w: mp.mpf(sp.N(v[i].subs({Delta: sp.Float(mp.nstr(D, 45), 45), J: Jv}), 45)) for i, w in enumerate(ws)}

X = [mp.mpf('1.3'), mp.mpf('0.4'), mp.mpf('-0.5')]; Y = [mp.mpf('0.9'), mp.mpf('-0.3'), mp.mpf('0.2')]
pr("# heavy vertex of <JJJJ>: max_h |exact/closed - 1| over the 16 helicity components, per (a,b); leakage; m*(sym-1)")
for (D0, J0) in (('50.25', 20), ('40.5', 8), ('60.5', 40)):
    for scale in (10, 100, 1000, 10000):
        D = mp.mpf(D0)*scale; Jv = J0*scale
        Ap, Am, m, beta = weyl(D, Jv)
        row = []
        R = {}
        for a in range(2):
            for b in range(2):
                cL = coeffs(vecs[a], D, Jv); cR = coeffs(vecs[b], 3 - D, -1 - Jv)
                res, leak, _ = helicity_components(1, cL, cR, D, Jv, X, Y, m)
                dev = 0
                for h, v in res.items():
                    pred = V_l1(a+1, Ap, Am, m, h[0], h[1])*V_l1(b+1, Ap, Am, m, h[2], h[3])
                    R[(a, b, h)] = v/pred
                    dev = max(dev, abs(v/pred - 1))
                row.append("(%d,%d): %s leak %s" % (a+1, b+1, mp.nstr(dev, 3), mp.nstr(abs(leak)/abs(pred), 2)))
        # Weyl antisymmetry of the O(1/m) term: R(a,b,(h1,h2,h3,h4)) + R(b,a,(h3,h4,h1,h2)) - 2 = O(1/m^2)
        sym = max(abs(R[(a, b, h)] + R[(b, a, (h[2], h[3], h[0], h[1]))] - 2) for (a, b, h) in R)
        pr(f"  Delta={mp.nstr(D, 8)} J={Jv} m={mp.nstr(m, 6)} beta={mp.nstr(beta, 4)}: " + "; ".join(row) + f";  m^2*max|sym-1| = {mp.nstr(m*m*sym/2, 4)}")
pr("# kernel conventions: K_{mu,0}(xi) vs F^nu_mu(xi) (e2e), and K_{0,mu}(xi)")
for (nu, xi) in ((0.0, 0.4), (1.0, 1.0), (3.0, 1.0)):
    for mu in (0, 2):
        k1 = kernel_mp(nu, mu, 0, xi); k2 = kernel_mp(nu, 0, mu, xi); f = Fnu_mp(nu, mu, xi)
        kn = kernel_np(nu, [(mu, 0)], [xi])[(mu, 0)][0]
        pr(f"  nu={nu} xi={xi} mu={mu}: K_mu0 = {mp.nstr(k1, 12)}  F = {mp.nstr(f, 12)}  K_0mu = {mp.nstr(k2, 12)}  numpy K_mu0 = {kn:.12g}")
open('out_t3_vertex_l1.txt', 'w').write("\n".join(out) + "\n")
