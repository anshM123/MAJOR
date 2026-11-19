"""v2 (H7_I14_LIGHTSIDE.md §13): (HN) -- model size of a spinning magic-block mismatch.

If the leading-Regge kernel of the spinning CLS-D functional were eta x (S (x) S') P_l with derivative-free transverse structures S, S'
(products of Lemma-L D1, D2 at each pair), then (exactly, since eta-multiplication commutes with S and Lemma L maps S P_l' -> M_S(nu') P_l'):
   preserving heavy image / [C^2 gamma^2 M_S M_S'](nu)  =  (1/2)[g rho P_{l-1} + g~ conj(rho) P_{l+1}],
   rho(nu) = [M_S M_S'](nu + i) / [M_S M_S'](nu),
i.e. the scalar magic-block coefficient g is replaced by g rho (conjugation-symmetric, so Re does not remove it); (HN) with the
scalar de-shift needs rho == 1.  This script tabulates |rho - 1| at a = 7 for S = S' in {D1^2, D1D2, D2^2,
D1, D2}.  Output: cache/v2_hn_model.txt
"""
import os, io
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
buf = io.StringIO()
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")
a = 7.0
def dd(nu, b): return nu*nu + (a - b)**2
def M(s, nu):
    d5, d9 = dd(nu, 2.5), dd(nu, 4.5)
    N1 = 64*a**4 - 768*a**3 + 64*a*a*nu*nu + 3280*a*a - 320*a*nu*nu - 5840*a + 16*nu**4 + 488*nu*nu + 3705
    N2 = 64*a**3 - 32*a*a*nu*nu - 584*a*a + 224*a*nu*nu + 1720*a - 16*nu**4 - 424*nu*nu - 1641
    N3 = 128*a*a - 128*a*nu*nu - 800*a + 16*nu**4 + 360*nu*nu + 1113
    return {'D1^2': N1/(64*d5*d9), 'D1D2': -N2/(64*d5*d9), 'D2^2': N3/(64*d5*d9),
            'D1': (nu*nu + 2*(a - 2)**2 + 0.25)/(2*d5), 'D2': (nu*nu - 2*a + 17/4)/(2*d5)}[s]
nus = (0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 25.0, 50.0)
pr("|rho(nu) - 1|, rho = [M_S(nu+i)/M_S(nu)]^2  (S = S' at both pairs), a = 7")
pr("   S      " + "".join(f"{v:>10}" for v in nus))
for s in ('D1', 'D2', 'D1^2', 'D1D2', 'D2^2'):
    r = [abs((M(s, v + 1j)/M(s, v))**2 - 1) for v in nus]
    pr(f"   {s:6s} " + "".join(f"{x:10.2e}" for x in r))
pr("arg rho / pi:")
for s in ('D1^2', 'D2^2'):
    r = [np.angle((M(s, v + 1j)/M(s, v))**2)/np.pi for v in nus]
    pr(f"   {s:6s} " + "".join(f"{x:+10.2e}" for x in r))
pr("large-nu: rho - 1 ~ 2i M_S'(nu)/M_S(nu) = O(nu^-3) for these structures (M_S -> const + O(nu^-2)).")
open(os.path.join(HERE, 'cache', 'v2_hn_model.txt'), 'w').write(buf.getvalue())
