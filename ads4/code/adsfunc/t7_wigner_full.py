"""t7: flat positivity of the graviton functionals with the FULL Wigner function d^J_{44} = (1 - p^2/m^2)^4 dtilde^J_{44} in the +- channel
(ADS4_FUNCTIONAL.md section 7: the density that the AdS moments carry if SKELETON_TEST.md section 8 (currents, reading B) carries over to
gravitons), versus the helicity-stripped dtilde^J_{44} of arXiv:2201.06602 (2.16).  ++ is unchanged (d_00 = dtilde_00).

Per +- state (normal form, u = 1/m^2, M = 1):
   original:  FF  = int [psi2 (2 - u p^2) - u psi3] dtilde44(1 - 2up^2) dp + a (2JJ - 41 - 1) u^3
   modified:  FF' = int [psi2 (2 - u p^2) - u psi3] (1 - u p^2)^4 dtilde44(1 - 2up^2) dp + a (2JJ - 33 - 1) u^3
(the forward constant: m^10 d_{p^2}C_4|_0 = D'(0) + 2m^2 D''... gives 41 -> 33 because d_{p^2}(1-p^2/m^2)^4|_0 = -4/m^2).
Functionals: F_g3 completed (a = 5) and literal (no B5: FW = a(2JJ - 41)u^3 -> a(2JJ - 33)u^3), M01, M02 (ads4/code/robust/functionals_certified.json).
Grid: m in [1,2] step 0.005, [2,40] step 0.1, [40,120] step 2; J = 4 ... min(200 m, 30000) (bhat <= 400).  Floating point; not a certificate.
"""
import json, os
import numpy as np
from fractions import Fraction as Fr
from numpy.polynomial.legendre import leggauss
from common import jac_all, PSI2, PSI3, peval

HERE = os.path.dirname(os.path.abspath(__file__))
fun = {'Fg3 completed': (PSI2, PSI3, 5.0, True), 'Fg3 literal': (PSI2, PSI3, 5.0, False)}
d = json.load(open(os.path.join(HERE, '..', 'robust', 'functionals_certified.json')))
for k in ['M01', 'M02']:
    fun[k] = (np.array([float(Fr(x)) for x in d[k]['psi2_exact']]), np.array([float(Fr(x)) for x in d[k]['psi3_exact']]),
              float(Fr(d[k]['a_exact'])), True)

x, w = leggauss(1600)
p = (x + 1) / 2
w = w / 2
MS = np.concatenate([np.arange(1.0, 2.0, 0.005), np.arange(2.0, 40.0, 0.1), np.arange(40.0, 120.01, 2.0)])

res = {k: {'orig': (np.inf, None), 'mod': (np.inf, None), 'orig_rel': (np.inf, None), 'mod_rel': (np.inf, None)} for k in fun}
for m in MS:
    u = 1 / m ** 2
    Jmax = int(min(30000, max(60, 200 * m)))
    D = jac_all(Jmax - 4, 0, 8, 1 - 2 * u * p ** 2)          # dtilde44 for J = 4 .. Jmax
    Js = np.arange(4, Jmax + 1)
    JJ = Js * (Js + 1.0)
    fac = (1 - u * p ** 2) ** 4
    for k, (c2, c3, a, comp) in fun.items():
        Psi = peval(c2, p) * (2 - u * p ** 2) - u * peval(c3, p)
        for tag, wt, c4 in [('orig', w * Psi, 41), ('mod', w * Psi * fac, 33)]:
            B = D @ wt
            FW = a * (2 * JJ - c4 - (1 if comp else 0)) * u ** 3
            F = B + FW
            N = np.abs(D) @ np.abs(wt) + np.abs(FW)
            i = np.argmin(F)
            if F[i] < res[k][tag][0]:
                res[k][tag] = (F[i], (round(m, 3), int(Js[i])))
            r = F / N
            i = np.argmin(r)
            if r[i] < res[k][tag + '_rel'][0]:
                res[k][tag + '_rel'] = (r[i], (round(m, 3), int(Js[i])))
for k in fun:
    print('%-14s +- channel:  original  min FF = %+.4e at (m,J)=%s, min FF/N = %+.3e at %s' % (k, *res[k]['orig'], *res[k]['orig_rel']))
    print('%-14s              full d44  min FF = %+.4e at (m,J)=%s, min FF/N = %+.3e at %s' % ('', *res[k]['mod'], *res[k]['mod_rel']))
# exact values at the threshold J = 4, m = 1 (dtilde44 = 1 for J = 4)
for k, (c2, c3, a, comp) in fun.items():
    Psi = peval(c2, p) * (2 - p ** 2) - peval(c3, p)
    o = np.sum(w * Psi) + a * (40 - 41 - (1 if comp else 0))
    mo = np.sum(w * Psi * (1 - p ** 2) ** 4) + a * (40 - 33 - (1 if comp else 0))
    print('%-14s J=4, m=1: original %.6f  full-d44 %.6f' % (k, o, mo))
