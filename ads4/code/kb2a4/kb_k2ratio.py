"""kb_k2ratio.py -- the k = 2 (EE, j'=0, f = v'-u') T-block light values, Einstein and |g3|^2, in the x0 = inf momentum-space reduction
(ads4/K_B2A4.md Sec. 4).

Rest-frame integrand (P = (1,0,0), mostly-minus; e_i = 1 - cos phi_i0, e13 = 1 - cos phi13):
   Psi^+  propto  int dphi dphi0  C0 (2/e13) { T31 d3d1F/F + (d3F/F)[T31 lam1 + tau1/2] + (d1F/F)[T31 lam3 + tau3/2] },
   C0 = (2e13)^{1+delta/2} (2e1)^{-delta/2} (2e3)^{-delta/2},  F = (z0.P)^alpha (P^2)^kappa, alpha = delta-1, kappa = -11/4 - i nu/2,
   dF/F = alpha e + 2 kappa ;  d3d1F/F = (alpha e1 + 2k)(alpha e3 + 2k) - alpha e1 e3 + 2 kappa (e13 - 2) ;  lam = 3/2 ;
   Einstein: T31 = 2, tau = 6 ;  W x W: T31 = cos(4phi)/8, tau = 3cos(4phi)/8  (kb_onepoint.py).
Hence B_W = (cos 4phi/16) B_E pointwise, B_E = 2Y, Y = alpha(alpha-1) e1e3 + (2 alpha kappa + 3 alpha)(e1+e3) + 2 kappa e13 + 4kappa^2 + 8kappa.
Every monomial integral is a three-point circle integral = Dixon's identity (closed Gamma form); cos 4phi = T4(1 - e13).
Outputs: e_EE(nu)/I0 (compared with the Lemma-R Einstein value of the EE functional, R_-D(nu+i)+R_+D(nu-i)), Theta(nu) = Psi[WW]/(t4^2 Psi[EE])
exactly, its top coefficient (flat check 1/(8*240^2)), and an independent numerical evaluation of Theta by Fourier sums (no Dixon).
Usage: python3 -B kb_k2ratio.py -> cache/kb_k2ratio.txt
"""
import os
import sympy as sp
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
d, nu, x = sp.symbols('delta nu x')


def rf(a, n):
    return sp.Mul(*[a + k for k in range(n)]) if n > 0 else sp.Integer(1)


def dixon_ratio(A0, B0, C0, a, b, c):
    """D(A0+c, B0+a, C0+b)/D(A0,B0,C0), D(A,B,C) = G(1+A+B+C)G(1+2A)G(1+2B)G(1+2C)/[G(1+A)G(1+B)G(1+C)G(1+A+B)G(1+B+C)G(1+C+A)]
    = (2pi)^-3 int |x1-x3|^{2A}|x1-x0|^{2B}|x3-x0|^{2C} on the unit circle (Dixon)."""
    r = rf(1+A0+B0+C0, a+b+c)*rf(1+2*A0, 2*c)*rf(1+2*B0, 2*a)*rf(1+2*C0, 2*b)
    return r/(rf(1+A0, c)*rf(1+B0, a)*rf(1+C0, b)*rf(1+A0+B0, c+a)*rf(1+B0+C0, a+b)*rf(1+C0+A0, b+c))


A0, B0, C0 = d/2, -d/2, -d/2
I = lambda a, b, c: sp.Rational(1, 2**(a+b+c))*dixon_ratio(A0, B0, C0, a, b, c)
al = d - 1
ka = -sp.Rational(11, 4) - (d - sp.Rational(1, 2))/2
Y = [(al*(al-1), 1, 1, 0), (2*al*ka+3*al, 1, 0, 0), (2*al*ka+3*al, 0, 1, 0), (2*ka, 0, 0, 1), (4*ka**2+8*ka, 0, 0, 0)]
e = sp.Symbol('e')
T4 = sp.Poly(sp.expand(8*(1-e)**4 - 8*(1-e)**2 + 1), e)


def main():
    out = []

    def pr(s):
        print(s, flush=True); out.append(s)
    sub = {d: sp.Rational(1, 2) + sp.I*nu}
    EE = sp.factor(sp.simplify(sp.together(sum(c*I(a, b, cc) for c, a, b, cc in Y)).subs(sub)))
    WW = sp.together(sum(c*v*I(a, b, cc+k[0]) for c, a, b, cc in Y for k, v in zip(T4.monoms(), T4.coeffs()))/16)
    Th = sp.factor(sp.simplify((WW.subs(sub))/(EE)))
    pr('# kb_k2ratio')
    pr('Einstein T-block, k=2 EE, rest-frame reduction:  Psi[EE] / (U(nu) I0(nu)) = %s' % EE)
    # Lemma-R Einstein value of the EE functional (H7 Sec. 5, EE leading-Regge coefficients R_-+ = R^sc_-+ rho(nu +- i))
    Rms = -(2*nu+sp.I)*(2*nu+11*sp.I)**2/(16*nu); Rps = Rms.subs(nu, -nu)
    rho = lambda z: 32*(4*z**2+1)*(4*z**2+9)/((4*z**2+25)*(4*z**2+81))
    D = lambda z: 1/(z**2 + sp.Rational(9, 4))
    ELR = sp.factor(sp.simplify(Rms*rho(nu+sp.I)*D(nu+sp.I) + Rps*rho(nu-sp.I)*D(nu-sp.I)))
    pr('Lemma-R Einstein value  R_- D(nu+i) + R_+ D(nu-i) = %s' % ELR)
    ratio = sp.factor(sp.simplify(ELR/EE))
    pr('ratio (must be free of real zeros/poles; Gamma-type normalisation U(nu) I0(nu)/N_com absorbs it): %s' % ratio)
    pr('Theta(nu) = Psi[WxW]/(t4^2 Psi[ExE]) = %s' % Th)
    ser = sp.series(Th.subs(nu, 1/x), x, 0, 2)
    pr('large-nu: Theta = %s' % ser)
    top = sp.limit(Th/nu**8, nu, sp.oo)
    pr('top coefficient = %s ; flat value |g3|^2 nu^8/8 with t4^2 = 240^2 |g3|^2 requires 1/(8*240^2) = %s : %s' % (
        top, sp.Rational(1, 8*240**2), top == sp.Rational(1, 8*240**2)))
    # h1-reading light function: per-nu Einstein (forward-normalised) value E_EE/(R_-+R_+), times 8*57600*Theta, over nu^6
    S = sp.factor(sp.simplify(Rms*rho(nu+sp.I) + Rps*rho(nu-sp.I)))
    g2 = sp.factor(sp.simplify(8*57600*Th*ELR/S/nu**6))
    pr('R_- + R_+ = %s' % S)
    pr('g2(nu) := [|g3|^2 light function, h1 reading]/flat = %s' % g2)
    pr('   large nu: g2 = %s   => a4 = -581/4' % sp.series(g2.subs(nu, 1/x), x, 0, 7))
    P = sp.factor(sp.simplify(g2*S*nu**6))
    pr('   (R_-+R_+) nu^6 g2 = P(nu) = %s  (regular on the real axis)' % P)
    # independent numerical check of Theta by Fourier sums (no Dixon, no T4 expansion)
    mp.mp.dps = 25

    def b(s, n):
        return mp.e**(-1j*mp.pi*n)*mp.gamma(2*s+1)*mp.rgamma(1+s-n)*mp.rgamma(1+s+n)

    def Lam(nv, s13, s10, s30, m, N=2500):
        return sum(b(s10, n)*b(s30, n)*b(s13, n - m) for n in range(-N, N))
    for nv in (mp.mpf('1.3'), mp.mpf(5)):
        dl = mp.mpf(1)/2 + 1j*nv
        aln = dl - 1; kan = -mp.mpf(11)/4 - 1j*nv/2
        Yn = [(aln*(aln-1), 1, 1, 0), (2*aln*kan+3*aln, 1, 0, 0), (2*aln*kan+3*aln, 0, 1, 0), (2*kan, 0, 0, 1), (4*kan**2+8*kan, 0, 0, 0)]
        num = 0; den = 0
        for c, a, bb, cc in Yn:
            s13, s10, s30 = dl/2 + cc, -dl/2 + a, -dl/2 + bb
            w = c*mp.mpf(2)**(-(a+bb+cc))
            den += w*Lam(nv, s13, s10, s30, 0)
            num += w*(Lam(nv, s13, s10, s30, 4) + Lam(nv, s13, s10, s30, -4))/2
        thn = num/den/16
        tha = complex(Th.subs(nu, sp.Float(str(nv), 30)).evalf(25))
        pr('Fourier-sum check nu=%s: Theta_num = %s, Theta_exact = %s, rel %.1e' % (nv, mp.nstr(thn, 14), mp.nstr(mp.mpc(tha), 14),
                                                                                  float(abs(thn/mp.mpc(tha) - 1))))
    return out


if __name__ == '__main__':
    o = main()
    open(os.path.join(HERE, 'cache', 'kb_k2ratio.txt'), 'w').write('\n'.join(o) + '\n')
