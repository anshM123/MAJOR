"""tailaudit/ta_cutoff_sup.py [certified, Arb]: the sup over [1, Lam+80] of the per-cell enclosure of |f_theta - f| <= 2|theta(nu-i)-1||g|
(same cells and enclosures as tailcert/tc_cutoff.py, which outputs only integrals and the sup beyond Lam+20); TAIL_CERT §4 quotes
sup|f_theta - f| <= 1e-21.  Also records where the sup sits.  Output: cache/ta_cutoff_sup.txt"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); ADS4 = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(ADS4, 'code', 'contour')); sys.dont_write_bytecode = True
import cb_common as K
from ew_common import horner, QRN, QRD, RMN, RMD
from flint import arb, acb, ctx
ctx.prec = 128
I = acb(0, 1)
W = K.load(); Lam = float(W.Lam.mid()); nuc = W.Lam + 100
h = 0.025; rho = 0.3; m = 1 + h; best = 0.0; where = None
while m < Lam + 80:
    z = acb(arb(m, rho), arb(0, rho)); zeta = z - I; re_hi = m + rho
    if re_hi <= float((nuc - 10).lower()):
        imz = arb(1 + rho)
        th1 = (imz*imz/100).exp()/2*(((arb(m - rho) + nuc)/10).erfc() + ((nuc - arb(re_hi))/10).erfc())
    else:
        th1 = 1 + (arb(1 + rho)**2/100).exp()
    n2 = zeta*zeta
    Y = (W.y[0] + W.y[1]*n2 + W.y[2]*n2*n2)*((-n2/2).exp())
    gsup = abs(horner(QRN, zeta)/horner(QRD, zeta))*W.B_disc_bound(m, rho) + abs(Y*(n2 + 1)*(n2 + 1).log()*horner(RMN, zeta)/horner(RMD, zeta))
    v = float((2*th1*gsup).upper())
    if v > best:
        best, where = v, m
    m += 2*h
s = f'max over cells in [1, Lam+80] of the enclosure of |f_theta - f| = {best:.3e} (cell centre nu = {where:.3f}); TAIL_CERT quotes <= 1e-21'
print(s); open(os.path.join(HERE, 'cache', 'ta_cutoff_sup.txt'), 'w').write(s + '\n')
