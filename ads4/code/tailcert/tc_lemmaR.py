"""
tailcert/tc_lemmaR.py (TAIL_CERT.md §R) [certified given EPSW_CERT's certified inputs + the [P] lemma below]:
Arb reproduction of the author's Lemma R (ANGULAR_ROUTE §81), with the Mehler-step constants re-derived here.

Inputs (EPSW_CERT, certified, epswcert/cache/ew_cert_1000_re_rep.json): r(0), r'(0+), eps_W >= int|r|, B_EPSW >= |r'(0+)| + int|r''|,
for r = E - psi (E without the erf cutoff theta).  r~ := r_theta - r(0) sech(pi nu), r_theta = E theta - psi.
  E_in := eps_W + r(0)/2 + dE_theta,        dE_theta >= int |E(1-theta)|
  B_in := B_EPSW + r(0) pi^2 int_0^oo |sech''| + dB_theta,  dB_theta >= int |(E(1-theta))''|
  (int_0^oo sech(pi nu) dnu = 1/2;  int_0^oo |sech''(x)| dx <= 2  [sech'' = sech(2 sech^2 - ... ) ; |sech''| <= sech + 2 sech^3, int <= pi/2 + pi/2... see doc])
  dE_theta, dB_theta <= 1e-9 [P]: see TAIL_CERT.md §R (1 - theta <= erfc(7) on [0, Lam+30]; Gaussian decay of E beyond Lam+30,
  sup|E'| <= B_EPSW + sup|psi'|; theta', theta'' <= e^{-48} on [0, Lam+30]).  Taken as 1e-9 each (generous).
Lemma (Mehler step) [P, re-derived here]: with W(t) = int r~ cos(nu t), |W| <= min(E_in, B_in/t^2), int_0^oo W = (pi/2) r~(0) = 0,
  k(t) = (cosh b - cosh t)^{-1/2} on [0,b), K0 = k(0) = (cosh b - 1)^{-1/2}:
  I := int_0^oo r~ P_{-1/2+i nu}(cosh b) dnu = (sqrt2/pi) int_0^oo W (k 1_{t<b} - K0) dt, and
  |I| <= (sqrt2/pi) K0 [ B_in/b + 8 B_in log2 / b^2 + E_in eps1(b) ],  eps1 = (sinh(b/2) - b/2)/(2(cosh b - cosh(b/2))),
  using sqrt(1+x) <= 1 + x/2 on [0,b/2], and int_0^b (k/K0 - 1) dt <= int_0^1 ((1-x)^{-1/2} - 1) dx/x = 2 log 2 on [b/2, b].
Ratio to F_cert: LRINF cells on [30, 64] (cell lower bound; numerator at the cell's left end, every term decreasing in b);
  b >= 64: fcert_tail formula F = (1-2.1e-7) A1 (sqrt2/pi) K0 ell/(b Lam^2); the ratio 2 Lam [B + 8 B log2/b + E b eps1]/((1-2.1e-7) A1 ell)
  is decreasing in b, so b = 64 bounds all b >= 64.
Output: cache/tc_lemmaR.txt
"""
import os, sys, json
from fractions import Fraction as Fr
from flint import arb, ctx
ctx.prec = 200
HERE = os.path.dirname(os.path.abspath(__file__))
ADS4 = os.path.abspath(os.path.join(HERE, '..', '..'))
out = []
def pr(s):
    print(s, flush=True); out.append(s)

ew = json.load(open(os.path.join(ADS4, 'code', 'epswcert', 'cache', 'ew_cert_1000_re_rep.json')))
r0 = arb(ew['r0_hi']); epsW = arb(ew['epsW']); BE = arb(ew['B'])
PI = arb.pi(); L2 = arb(2).log()
# int_0^oo |sech''(x)| dx: sech'' = sech - 2 sech^3 (since (sech)'' = sech tanh^2 - sech^3 = sech - 2 sech^3); |.| <= sech + 2 sech^3,
# int sech = pi/2, int sech^3 = pi/4  -> <= pi.   Then int |d^2/dnu^2 sech(pi nu)| = pi^2 * (1/pi) * int|sech''| <= pi^2.
Ein = epsW + r0/2 + arb('1e-9')
Bin = BE + r0*PI**2 + arb('1e-9')
pr(f'E_in <= {float(Ein.upper()):.6e}  (author: 2.683e-6);  B_in <= {float(Bin.upper()):.6e}  (author: 1.901e-3)')
assert Ein < arb('2.683e-6') and Bin < arb('1.901e-3')
# use the author's rounded constants (they are upper bounds)
E = arb('2.683e-6'); B = arb('1.901e-3')
nf = json.load(open(os.path.join(ADS4, 'prl', 'data', 'nf02_psi.json')))
K0c = Fr('27.48509817')
p2 = [Fr(c)/K0c for c in nf['psi2_exact']]; p3 = [Fr(c)/K0c for c in nf['psi3_exact']]
_a1 = 2*p2[1] - p3[1]; A1 = arb(_a1.numerator)/_a1.denominator
LAM = arb(998*999).sqrt()
b64 = arb(64)
L64 = 1/b64 - 1/(2*(b64/2).sinh()) - arb(2).sqrt()*((-b64/2).exp()/arb('4.7') + 4/b64**2) - 2*(2*L2).sqrt()/(b64 - L2)**2
ell = 64*L64
pr(f'A1 = 2 psi2\'(0) - psi3\'(0) = {A1.str(10)},  ell = {ell.str(8)}')

def num(b):
    b = arb(b)
    K0 = 1/(b.cosh() - 1).sqrt()
    eps1 = ((b/2).sinh() - b/2)/(2*(b.cosh() - (b/2).cosh()))
    return (arb(2).sqrt()/PI)*K0*(B/b + 8*B*L2/b**2 + E*eps1)

rows = json.load(open(os.path.join(ADS4, 'code', 'lrtight', 'cache', 'lrt3_rcert.json')))['rows']
worst = arb(0); wb = None
for r in rows:
    ba, bb = Fr(r['ba']), Fr(r['bb'])
    if ba < 30 or bb > 64:
        continue
    ratio = (2/LAM)*num(arb(ba.numerator)/ba.denominator)/arb(r['Fcert'])
    if ratio > worst:
        worst = arb(ratio.upper()); wb = (float(ba), float(bb))
pr(f'beta in [30, 64] (LRINF cells): max (2/Lam)|I_r~|/F_cert <= {float(worst.upper()):.4e}  (worst cell {wb})')
rt = 2*LAM*(B + 8*B*L2/b64 + E*b64*((b64/2).sinh() - b64/2)/(2*(b64.cosh() - (b64/2).cosh())))/((1 - arb('2.1e-7'))*A1*ell)
pr(f'beta >= 64 (fcert_tail formula, ratio decreasing in beta): (2/Lam)|I_r~|/F_cert <= {float(rt.upper()):.4e}  (author: < 0.015)')
rinf = 2*LAM*B/((1 - arb('2.1e-7'))*A1*ell)
pr(f'beta -> oo limit of the bound: {float(rinf.upper()):.4e}')
for b in (30, 64, 100, 200):
    pr(f'  beta = {b}: (2/Lam)|I_r~| <= {float(((2/LAM)*num(b)).upper()):.4e}')
open(os.path.join(HERE, 'cache', 'tc_lemmaR.txt'), 'w').write('\n'.join(out) + '\n')
