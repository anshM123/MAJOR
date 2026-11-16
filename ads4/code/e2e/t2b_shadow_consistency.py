"""
Consistency of the shadow coefficients b^{(Delta_34)} with conservation.
The right structure of the spinning block in the Lorentzian shadow representation is
     R~ = sum_w c_w(Delta,J) / b^{(D34_w)}_{Delta,J} * D_w S3^S(Delta_3 + #a, Delta_4 + #b),      D34_w = (l+1+#a) - (Dphi+#b),
(S3^S: scalar structure with (Delta,J) -> (3-Delta, -1-J)).  The shadow transform preserves conservation and the conserved
parity-even structure is unique, hence   c_w(Delta,J)/b^{(D34_w)} = kappa * c_w(3-Delta,-1-J)   for all words w.
Check: kappa_w := c_w(Delta,J) / (b_w c_w(3-Delta,-1-J)) must be independent of w.  This tests the Delta_34-dependence of b
(the Gamma-function generalisation of CHMRS (3.6)) exactly, and gives kappa (the shadow coefficient of the conserved structure).
"""
import mpmath as mp, sympy as sp
from structures import get as get_struct, Dphi as symDphi, WORDS
from embed import Delta as symDelta, J as symJ
from formulas import log_b
mp.mp.dps = 30
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
def kappas(l, D, Jv, Dp):
    cs, _ = get_struct(l)
    res = {}
    for w in WORDS[l]:
        c1 = cs[w].subs({symDelta: sp.nsimplify(D), symJ: Jv, symDphi: sp.nsimplify(Dp)})
        c2 = cs[w].subs({symDelta: 3 - sp.nsimplify(D), symJ: -1 - Jv, symDphi: sp.nsimplify(Dp)})
        D34 = (l + 1 + w.count('a')) - (Dp + w.count('b'))
        res[w] = mp.mpf(sp.Rational(c1, c2).p)/mp.mpf(sp.Rational(c1, c2).q)*mp.exp(-log_b(mp.mpf(D), Jv, D34))
    return res
pr("# kappa_w * b^{(D34 of the last word)}  for each word (must coincide);  last column: the same computed with b^{(0)} for all words (must NOT coincide)")
for l in (1, 2):
    for (D, Jv, Dp) in ((9.75, 3, 1.75), (30.25, 6, 1.75), (30.25, 6, 2.5), (100.5, 40, 1.75)):
        k = kappas(l, D, Jv, Dp)
        wl = WORDS[l][-1]
        D34l = (l + 1 + wl.count('a')) - (Dp + wl.count('b'))
        bl = mp.exp(log_b(mp.mpf(D), Jv, D34l))
        row = [mp.nstr(k[w]*bl, 20) for w in WORDS[l]]
        cs, _ = get_struct(l)
        row0 = []
        for w in WORDS[l]:
            sub1 = {symDelta: sp.nsimplify(D), symJ: Jv, symDphi: sp.nsimplify(Dp)}
            sub2 = {symDelta: 3 - sp.nsimplify(D), symJ: -1 - Jv, symDphi: sp.nsimplify(Dp)}
            row0.append(mp.nstr(mp.mpf(str(sp.N(cs[w].subs(sub1)/cs[w].subs(sub2), 30))), 8))
        pr(f"  l={l} Delta={D} J={Jv} Dphi={Dp}: {row}   | b^(0): {row0}")
open('out_t2b_shadow_consistency.txt', 'w').write("\n".join(out) + "\n")
