"""
t1: structure of the exact (CHMRS order-8) heavy moments in the bulk-point regime and the size of the genuinely-AdS remainder.

(a) seriesCfixedJ = m^{2k} C_{k,nu}[Delta,J] at fixed J, d = 3:  2 + sum_j Q_{k,2j}(nu, J2)/m^{2j}.  Split each Q_{k,2j} by degree in nu^2:
    degree j   = leading-nu part (must equal the expansion of the bulk-point formula; checked symbolically);
    degree j-1 = next-to-leading: R1_j(J2) nu^{2j-2} -> the 1/m^2 remainder function  r_1(x) = sum_j R1_j x^{2j-2}  (x = nu/m);
    degree j-2 = NNLO: r_2(x) = sum_j R2_j x^{2j-4} (1/m^4 remainder).
    So  m^{2k} C = BP(x) + r_1(x)/m^2 + r_2(x)/m^4 + ...  (r_1 known to x^6, r_2 to x^4).
(b) the same for the forward moments at nu = 0 (all orders available exactly): C_{4}, d/dnu^2 C_4, C_5 at nu = 0.
(c) checks: S_nu formula (CHMRS app. E) -> BP at small J and -> Regge leading at small nu (numerically).
Output: cache/out_t1.txt
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sympy as sp, numpy as np, mpmath as mp
from chmrs_series import load, d, k, Df, nu, J2
import skel

out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)

s = load()
x, eps = sp.symbols('x eps', positive=True)
cs = [sp.expand(c.subs(d, 3)) for c in s['seriesCfixedJ']]
pr("# t1 (a): seriesCfixedJ split by degree in nu^2  (d = 3)")
# bulk point expansion: (2-x^2)/(1-x^2)^{1+k/2} P_J(1-2x^2), with P_J(1-2x^2) = 2F1(-J, J+1; 1; x^2) = sum_n (-J)_n (J+1)_n/n!^2 x^{2n}
# (-J)_n (J+1)_n is a polynomial in J2 = J(J+1):  prod_{i=0}^{n-1} (i(i+1) - J2)
def PJ_series(N):
    return sum(sp.prod([sp.Integer(i*(i + 1)) - J2 for i in range(n)])/sp.factorial(n)**2*x**(2*n) for n in range(N + 1))
BP = sp.series((2 - x**2)*(1 - x**2)**(-1 - k/2), x, 0, 10).removeO()
BP = sp.expand(sp.series(sp.expand(BP*PJ_series(4)), x, 0, 10).removeO())
ok = True
for j in range(5):
    poly = sp.Poly(cs[j], nu)
    lead = sp.expand(poly.coeff_monomial(nu**(2*j)))
    bp = sp.expand(BP.coeff(x, 2*j))
    diff = sp.simplify(lead - bp)
    pr(f"  order 1/m^{2*j}: leading nu^{2*j} coefficient - bulk-point expansion = {diff}")
    ok = ok and diff == 0
pr(f"  => leading-nu terms of all Q_(k,2j), j <= 4, reproduce the bulk-point formula exactly: {ok}")

R1 = [sp.factor(sp.Poly(cs[j], nu).coeff_monomial(nu**(2*j - 2))) for j in range(1, 5)]
R2 = [sp.factor(sp.Poly(cs[j], nu).coeff_monomial(nu**(2*j - 4))) for j in range(2, 5)]
pr("# NLO coefficients R1_j (1/m^2 remainder r_1(x) = sum_j R1_j x^(2j-2)), as functions of k, Df, J2:")
for j, r in enumerate(R1, 1):
    pr(f"  R1_{j} = {r}")
pr("# NNLO coefficients R2_j (1/m^4 remainder r_2(x) = sum_j R2_j x^(2j-4)):")
for j, r in enumerate(R2, 2):
    pr(f"  R2_{j} = {sp.expand(r)}")
import pickle
pickle.dump({'R1': [sp.srepr(r) for r in R1], 'R2': [sp.srepr(r) for r in R2]}, open(os.path.join(HERE, 'cache', 'remainder_coeffs.pkl'), 'wb'))

pr("# numbers: Df = 3 (skeleton of <TTTT>: a = 2 Df + k - 1 = k + 5), relative to the flat value 2 at x = 0")
for kk in (2, 3, 4, 5):
    for JJ in (0, 2, 4, 8):
        r1 = [float(r.subs({k: kk, Df: 3, J2: JJ*(JJ + 1)})) for r in R1]
        r2 = [float(r.subs({k: kk, Df: 3, J2: JJ*(JJ + 1)})) for r in R2]
        pr(f"  k={kk} J={JJ:2d}: R1 = {[round(v, 3) for v in r1]}   R2 = {[round(v, 1) for v in r2]}   ratios R1_(j+1)/R1_j = {[round(r1[i+1]/r1[i], 3) if r1[i] else None for i in range(3)]}")
pr("# same for Df = 1.75 and Df = 2 (other skeletons), k = 2:")
for DD in (1.75, 2.0):
    for JJ in (0, 2, 4):
        r1 = [float(r.subs({k: 2, Df: DD, J2: JJ*(JJ + 1)})) for r in R1]
        pr(f"  Df={DD} k=2 J={JJ}: R1 = {[round(v, 3) for v in r1]}")

# (b) forward moments at nu = 0 from seriesC (fixed eta, all J): exact coefficients
pr("# t1 (b): forward moments at nu = 0 (seriesC, fixed eta_AdS), Df = 3")
cC = [sp.expand(c.subs(d, 3)) for c in s['seriesC']]
from chmrs_series import al, eta, P
pr("  structure: m^{2k}C_{k,0} = 2P(eta) + sum_j c_j(eta; P,P',..)/m^{2j};  coefficient growth (CHMRS): c_j ~ eta^{j+n} P^(n)")
for kk in (4, 5):
    for et in (1.0, 1.02, 1.2, 1.54, 2.0, 3.76):
        t = skel.series_terms('seriesC', kk, 3.0, 0.0, et)
        pr(f"  k={kk} eta={et:5.2f}: terms (m^0..m^-8 coefficients) = {[float('%.4g' % v) for v in t]}")

# (c) S_nu formula checks
pr("# t1 (c): S_nu (CHMRS E.2-E.3) vs bulk-point (fixed J, nu ~ m) and vs Regge leading (nu << m)")
for (Delta, J, nn) in ((200.0, 2, 100.0), (200.0, 2, 150.0), (800.0, 4, 400.0), (2000.0, 0, 1500.0), (400.0, 100, 2.0), (400.0, 100, 10.0), (4000.0, 1000, 20.0)):
    m2, eta_ = skel.m2_eta(Delta, J)
    for kk in (2, 4):
        S = skel.S_formula(kk, nn, Delta, J)
        bp = skel.C_BP(kk, nn, m2, J)
        rg = skel.C_regge_lead(nn, eta_)
        pr(f"  Delta={Delta:7.1f} J={J:4d} nu={nn:7.1f} x={nn/np.sqrt(m2):.3f} k={kk}:  S={S: .6e}  BP={bp: .6e}  Regge2P={rg: .6e}   S/BP={S/bp if bp else np.nan:.5f}  S/Regge={S/rg if rg else np.nan:.5f}")
open(os.path.join(HERE, 'cache', 'out_t1.txt'), 'w').write("\n".join(out) + "\n")
