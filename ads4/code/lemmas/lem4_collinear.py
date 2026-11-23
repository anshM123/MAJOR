r"""
lem4_collinear (LEMMAS_CLOSE.md §4): the a = b (collinear) set against the celestial kernel C_nu.
 (1) [P, sympy] At the joint light-ray endpoint t1, t3 -> 0 (a = 1/t1, b = 1/t3 -> inf) of the og_free gauge, the cross-ratio angle
     eta := (z + zb)/(2 sqrt(z zb)) satisfies eta^2 -> (a+b)^2/(2 y a b) >= 2/y = 1/sin^2(dphi/2): the collinear limit dphi -> 0 of the joint
     corner IS the large-eta (light-cone-edge) direction of the u-channel Regge box.
 (2) Uniform kernel bound |C_nu(dphi)| <= A13^{5/4} Kbar(dphi), Kbar ~ 2 log(1/dphi) + k0 (lem3 (a)); hence a light function
     |L| <= C dphi^{-p} pairs absolutely with C_nu, uniformly in nu, iff p < 7/2 (measure d(dphi) transverse to the diagonal).
     Table of the margin for the exponents met in the repo.
Output: cache/lem4_collinear.txt (seconds).
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import sympy as sp
import lem3_k2nu_lrc as L3

OUT = []


def pr(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); OUT.append(s)


if __name__ == '__main__':
    pr('# lem4_collinear: a = b integrability, LEMMAS_CLOSE.md §4')
    a, b, y, lam = sp.symbols('a b y lam', positive=True)
    S = 1 - (a + b)/(2*y); P = (2 - a)*(2 - b)/(8*y)
    eta2 = sp.simplify(S**2/(4*P))
    lim = sp.limit(eta2.subs(b, lam*a), a, sp.oo)
    pr('  (1) eta^2 = (z+zb)^2/(4 z zb) at a, b -> inf with b = lam a:', sp.simplify(lim), ';  minimum over lam at lam = 1:',
       sp.simplify(lim.subs(lam, 1)), ' (= 1/sin^2(dphi/2) with y = 2 sin^2(dphi/2))')
    for t, d in ((1e-3, 0.1), (1e-5, 0.01), (1e-7, 0.001)):
        aa = 1/t; yy = 2*np.sin(d/2)**2
        e = (1 - aa/yy)/(2*np.sqrt((2 - aa)**2/(8*yy) + 0j))
        pr(f'      t1 = t3 = {t:.0e}, dphi = {d}: |eta| = {abs(e):.6e}, 1/sin(dphi/2) = {1/np.sin(d/2):.6e}')
    pr('  (2) kernel bound and the integrability margin  (|C_nu| <= A13^{5/4} Kbar, A13 ~ dphi^2)')
    for d in (1e-4, 1e-3, 1e-2):
        kb = L3.K_trap(d, 0.0, power=-0.25).real
        pr(f'      dphi = {d:.0e}: Kbar = {kb:.4f}, 2 log(1/dphi) + 8.69 = {2*np.log(1/d) + 8.69:.4f};  A13^(5/4) Kbar = {(2 - 2*np.cos(d))**1.25*kb:.3e}')
    for lab, p in (('J0 < 4 tree terms (Lemma A: zero off a = b)', None), ('single free-scalar family (OG §3.6)', 3.0),
                   ('Witten-family local slopes (EIN §3.6; T-prime has Delta = 5, its own kernel differs: illustrative only)', 3.43), ('raw k = 4 Regge-edge growth eta^{5/2} (heuristic, §4)', 2.5)):
        if p is None:
            pr(f'      {lab}: integrand 0, nothing to integrate')
        else:
            pr(f'      {lab}: p = {p}: integrand dphi^({2.5 - p:+.2f}) x log -> {"integrable" if p < 3.5 else "NOT integrable"}, margin {3.5 - p:.2f}')
    with open(os.path.join(HERE, 'cache', 'lem4_collinear.txt'), 'w') as f:
        f.write('\n'.join(OUT) + '\n')
