"""
tailcert/tc_cutoff.py (TAIL_CERT.md §R.2, audit item ANGULAR_ROUTE §85): the erf-cutoff correction in the CORRECT form.
f_theta(nu) = 2 Re[theta(nu - i) g(nu)], g(nu) = chi(nu - i) R_-(nu - i);  f = 2 Re g (EPSW's E).  Hence
   |f_theta - f| <= 2 |theta(nu - i) - 1| |g(nu)|      (NOT via |Re g|),
   |(f_theta - f)''| <= 2 |((theta(. - i) - 1) g)''(nu)| <= 4 sup_{|z - m| <= rho} |(theta(z-i) - 1) g(z)| / (rho - h)^2   (Cauchy).
Certified on cells of width 2h = 0.05 over [1, Lam + 80] (discs of radius rho = 0.3 about the cell centres m):
   |theta(zeta) - 1| <= (1/2) e^{(Im zeta/10)^2} [erfc((Re zeta + nu_c)/10) + erfc((nu_c - Re zeta)/10)]  (Re zeta <= nu_c - 10;
      |erfc(x+iy)| <= e^{y^2} erfc(x), x > 0), else <= 1 + e^{(Im zeta)^2/100};
   |g(z)| <= |QR_-(zeta)| * (RBF disc bound, ew_common.Weight.B_disc_bound) + |Y l R_-(zeta)|  (Arb box enclosure, zeta = z - i).
Beyond Lam + 80: CONTOUR's Gaussian tail formula (cb_lines.tail, Y = 0) for int |g|, with |theta - 1| <= 2.02.
On [0, 1]: |theta(nu - i) - 1| <= e^{-12000}; the 1/nu parts of (g log)'' cancel in the real part (F2'(0) imaginary, theta(-i) real) [P],
so the [0,1] contributions are < e^{-11000} (not ball-evaluated).
Outputs: dE = int |f_theta - f|, dB = int |(f_theta - f)''|, dE_w = int e^{0.02 nu}|f_theta - f| (author's §85 quantity), sup_{nu>=Lam+20} 2|theta-1||g|.
"""
import os, sys, math, json
from flint import arb, acb, ctx
ctx.prec = 128
HERE = os.path.dirname(os.path.abspath(__file__))
ADS4 = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(ADS4, 'code', 'contour'))
import cb_common as K
import cb_lines as CL
from ew_common import horner, QRN, QRD, RMN, RMD
I = acb(0, 1)


def up(x):
    return float(arb(x).upper())


def main():
    W = K.load()
    Lam = float(W.Lam.mid()); nuc = W.Lam + 100
    h = 0.025; rho = 0.3
    dE = arb(0); dB = arb(0); dEw = arb(0); supfar = 0.0
    m = 1 + h
    end = Lam + 80
    while m < end:
        z = acb(arb(m, rho), arb(0, rho))          # box containing the disc |z - m| <= rho
        zeta = z - I
        re_hi = m + rho
        if re_hi <= float((nuc - 10).lower()):
            imz = arb(1 + rho)
            th1 = (imz*imz/100).exp()/2*(((arb(m - rho) + nuc)/10).erfc() + ((nuc - arb(re_hi))/10).erfc())
        else:
            th1 = 1 + (arb(1 + rho)**2/100).exp()
        QR = horner(QRN, zeta)/horner(QRD, zeta)
        n2 = zeta*zeta
        Y = (W.y[0] + W.y[1]*n2 + W.y[2]*n2*n2)*((-n2/2).exp())
        l = (n2 + 1)*(n2 + 1).log()
        Rm = horner(RMN, zeta)/horner(RMD, zeta)
        gsup = abs(QR)*W.B_disc_bound(m, rho) + abs(Y*l*Rm)
        hs = arb(th1*gsup)
        dE += 2*hs*(2*h)
        dEw += 2*hs*(2*h)*arb(0.02*(m + h)).exp()
        dB += 4*hs/(arb(rho) - h)**2*(2*h)
        if m >= Lam + 20:
            supfar = max(supfar, up(2*hs))
        m += 2*h
    intf, _ = CL.tail(W, 0.0, int(end) + 1)
    # CL.tail bounds int |f_theta| with the theta factor <= e^{1/100}; int |g| <= intf/2 (each of the two terms) times e^{0}
    tailg = intf
    dE += arb(2.02)*tailg; dEw += arb(2.02)*tailg*arb(0.02*(end + 1)).exp()*10   # (tail decays like e^{-(t-c)^2/1.62}: x10 covers e^{0.02 t})
    out = [f'dE  = int |f_theta - f|          <= {up(dE):.3e}',
           f'dEw = int e^(0.02 nu)|f_theta - f| <= {up(dEw):.3e}   (author §85: < 1e-7)',
           f'dB  = int |(f_theta - f)\'\'|       <= {up(dB):.3e}  (+ tail: second derivatives of the Gaussian tail, <= e^-300)',
           f'sup_(nu >= Lam+20) |f_theta - f|  <= {supfar:.3e}']
    for s in out:
        print(s, flush=True)
    open(os.path.join(HERE, 'cache', 'tc_cutoff.txt'), 'w').write('\n'.join(out) + '\n')
    json.dump(dict(dE=up(dE), dEw=up(dEw), dB=up(dB), supfar=supfar), open(os.path.join(HERE, 'cache', 'tc_cutoff.json'), 'w'))


if __name__ == '__main__':
    main()
