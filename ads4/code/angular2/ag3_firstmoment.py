"""ANGULAR_GO Sec. 4: scalar-exchange (J = 0) first-moment transition <O|E^(1)|T,+2> / <O|E|T,+2>, both parities, symbolic in Delta.
Uses code/kernels/pipeline.py + structures.py + kernels_lib.py READ-ONLY (KERNELS_CACHE redirected to angular2/cache, so nothing is
written in code/kernels).  First moment in the light-ray parameter alpha of pipeline.py (X13 = -alpha, X14 = -(alpha + 2 z.x)):
    I1^mid(a,b) = int (alpha - alpha_mid) X14^-a X13^-b = (a-b)/(2(a+b-2)) (2 z.x) I0(a,b),   alpha_mid = -z.x          [P]
per monomial: a = e14 - j, b = e13 - i, a - b = 2s-5+i-j, a+b-2 = esum-2-i-j; multiplying by (2 z.x) = 2 (z.x) is  l -> l+1, cf *= 2.
Variants reported: rule R1 = analytic (cf *= (a-b)/(a+b-2), l+1);  rule R2 = the coordinator's literal rule (cf *= (a-b)/(2(a+b-2)), l+1),
which is R1/2 (it equals R1 if the retarded time is u = alpha/2);  origin-centred moment  I1^0 = I1^mid + alpha_mid I0 = I1^mid - (z.x) I0.
Usage: python3 -B ag3_firstmoment.py  -> cache/ag3_firstmoment.txt"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ['KERNELS_CACHE'] = os.path.join(HERE, 'cache')
sys.path.insert(0, os.path.join(HERE, '..', 'kernels'))
import sympy as sp
import pipeline as pl          # read-only
import kernels_lib as kl       # read-only (cache redirected)
from ag_common import write

out = []
def pr(x):
    print(x, flush=True); out.append(str(x))

S = sp.Symbol('s'); Dl = sp.Symbol('Delta')

def shifted(expr, parity, mode):
    esum = 6 if parity == 'odd' else 5
    new = {}
    for (i, j, k, l, n1, n2, eps), cf in expr.items():
        if mode == 'origin':          # - (z.x) I0
            fac = pl.C(-1)
        else:
            num = 2*pl.s - 5 + i - j
            den = (esum - 2 - i - j) * (1 if mode == 'R1' else 2)
            fac = num*pl.C(1, den)
        key = (i, j, k, l + 1, n1, n2, eps)
        new[key] = new.get(key, pl.ZERO) + pl.redI(cf*fac)
    return {kk: v for kk, v in new.items() if v != 0}

for parity in ('even', 'odd'):
    monos, basis = kl.conserved(0, parity)
    pr('parity %s: %d monomials, %d conserved structure(s)' % (parity, len(monos), len(basis)))
    kin = pl.Kin('plus')
    e0 = [pl.scri_mono(m, kin) for m in monos]
    e1 = {md: [shifted(e, parity, md) for e in e0] for md in ('R1', 'R2', 'origin')}
    K = pl.Kbounds(e0 + sum(e1.values(), []), parity)
    for bi, coeffs in enumerate(basis):
        H = kl.poly_to_sym(kl.combine([pl.amplitude(e, kin, parity, K=K)[0] for e in e0], coeffs))
        pr('  structure %d: H (common prefactor G stripped) = %s' % (bi, sp.factor(H)))
        res = {}
        for md in ('R1', 'R2', 'origin'):
            M = kl.poly_to_sym(kl.combine([pl.amplitude(e, kin, parity, K=K)[0] for e in e1[md]], coeffs))
            r = sp.factor(sp.simplify(M/H))
            res[md] = r
            pr('    M_%s / H = %s   [Delta = 2s:  %s]' % (md, r, sp.factor(sp.simplify(r.subs(S, Dl/2)))))
        for md in ('R1', 'R2'):
            tot = sp.factor(sp.simplify((res[md] + res['origin']).subs(S, Dl/2)))
            pr('    origin-centred (%s + translation piece) / H = %s' % (md, tot))
        pr('    authors: M_mid/H = i(Delta^2-7)/6 ; check R1: %s ; R2: %s' % (
            sp.simplify(res['R1'].subs(S, Dl/2) - sp.I*(Dl**2 - 7)/6) == 0, sp.simplify(res['R2'].subs(S, Dl/2) - sp.I*(Dl**2 - 7)/6) == 0))
write('ag3_firstmoment.txt', out)
