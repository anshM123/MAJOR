"""
tailcert/tc_lownu_dv.py (TAIL_CERT.md §item3) [DV, mpmath]: SIGNED low-frequency residual contribution to the leading-Regge
(m = 0) term at large beta, relative to F_cert (coordinator request).
  G_r^low(beta) := (2/Lam) int_0^{nu1} r(nu) P_{-1/2+i nu}(cosh beta) dnu,   r = E - psi2(nu/Lam)/k0,  E(nu) = 2 Re[chi R_-](nu - i)
  (repaired weight k2_continuum_D_1000_rep.json, chi = chi_re; theta = 1 to 1e-200 here).  This is the low-nu piece of
  (1-c)^3 G^r(c) in FULLPOS3 §4.1(c) / FULLPOS5 §3 (whose absolute bound is (2/Lam) int|r| P_{-1/2}(cosh beta)).
  P_{-1/2+i nu}(cosh beta) = 2 Re[B(nu) w^{1/2-i nu} 2F1(1/2-i nu, 1/2-i nu; 1-2i nu; w)],  w = sech^2(beta/2),  B = Gamma(2i nu)/Gamma(1/2+i nu)^2.
  F_cert: LRINF cell values (beta < 64) and fp3_ht_large.fcert_tail formula (beta >= 64).
Output: cache/tc_lownu_dv.txt
"""
import os, sys, json, math
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
ADS4 = os.path.abspath(os.path.join(HERE, '..', '..'))
mp.mp.dps = 40
I = mp.mpc(0, 1)
d = json.load(open(os.path.join(ADS4, 'prl', 'data', 'k2', 'k2_continuum_D_1000_rep.json')))
CS = [(mp.mpf(c), mp.mpf(s)) for c, s in d['centres_sigma']]
X = [mp.mpf(float(v)) for v in d['x_re']]
Y = [mp.mpf(float(v)) for v in d['y_re']]
nf = json.load(open(os.path.join(ADS4, 'prl', 'data', 'nf02_psi.json')))
K0 = mp.mpf('27.48509817')
P2 = [mp.mpf(mp.fraction(*[int(t) for t in str(c).split('/')])) if '/' in str(c) else mp.mpf(c) for c in nf['psi2_exact']]
LAM = mp.sqrt(mp.mpf(998*999))
NEAR = [(c, s, x) for (c, s), x in zip(CS, X) if c < 40]


def Rm(nu):
    return -2*(2*nu-I)*(2*nu+I)**2*(2*nu+3*I)*(2*nu+5*I)*(2*nu+11*I)/(nu*(2*nu-7*I)*(2*nu-3*I)*(2*nu+7*I))


def chi(nu):
    Q = nu**2/(nu**2+mp.mpf(1)/4)**2
    Bg = 0
    for c, s, x in NEAR:
        g = lambda z: mp.exp(-z**2/(2*s*s))
        Bg += x*(g(nu-c) + (g(nu+c) if c > 0 else 0))
    n2 = nu*nu
    return Q*Bg + (Y[0] + Y[1]*n2 + Y[2]*n2*n2)*(n2+1)*mp.log(n2+1)*mp.exp(-n2/2)


def psi(nu):
    p = nu/LAM
    return sum(c*p**k for k, c in enumerate(P2))/K0


def r(nu):
    return 2*mp.re(chi(nu - I)*Rm(nu - I)) - psi(nu)


def Pcon(nu, beta):
    w = 1/mp.cosh(beta/2)**2
    Bn = mp.gamma(2*I*nu)/mp.gamma(mp.mpf(1)/2 + I*nu)**2
    a = mp.mpf(1)/2 - I*nu
    return 2*mp.re(Bn*w**a*mp.hyp2f1(a, a, 1 - 2*I*nu, w))


def fcert(beta):
    from fractions import Fraction as F
    rows = json.load(open(os.path.join(ADS4, 'code', 'lrtight', 'cache', 'lrt3_rcert.json')))['rows']
    if beta < 64:
        for rr in rows:
            if F(rr['ba']) <= beta < F(rr['bb']) or (beta == 64 and F(rr['bb']) == 64):
                return mp.mpf(rr['Fcert'])
    b = mp.mpf(beta); b64 = mp.mpf(64)
    L64 = 1/b64 - 1/(2*mp.sinh(b64/2)) - mp.sqrt(2)*(mp.exp(-b64/2)/mp.mpf(4.7) + 4/b64**2) - 2*mp.sqrt(2*mp.log(2))/(b64 - mp.log(2))**2
    p2 = [mp.mpf(c) for c in P2]
    p3 = [mp.mpf(c) for c in nf['psi3_exact']]
    A1 = (2*P2[1] - mp.mpf(nf['psi3_exact'][1]))/K0
    K0b = 1/mp.sqrt(mp.cosh(b) - 1)
    return (1 - mp.mpf(2.1e-7))*A1*(mp.sqrt(2)/mp.pi)*K0b*64*L64/(b*LAM**2)


if __name__ == '__main__':
    out = []
    def pr(s):
        print(s, flush=True); out.append(s)
    pr(f'r(0) = {mp.nstr(r(mp.mpf(0)), 6)}  (EPSW: +9.692e-14);  r(0.05)={mp.nstr(r(mp.mpf("0.05")),6)} r(0.31)={mp.nstr(r(mp.mpf("0.31")),6)}')
    b0 = mp.mpf(30)
    pr(f'check P at beta=30, nu=0.2: branch {mp.nstr(Pcon(mp.mpf("0.2"), b0), 15)}  legenp {mp.nstr(mp.legenp(-0.5+0.2j, 0, mp.cosh(b0), type=3).real, 15)}')
    for nu1 in (mp.mpf('0.31'), mp.mpf('1.0')):
        pts = [mp.mpf(0)] + [nu1*k/16 for k in range(1, 17)]
        rabs = mp.quad(lambda v: abs(r(v)), pts)
        pr(f'nu1 = {nu1}: int_0^nu1 |r| = {mp.nstr(rabs, 6)}')
        for beta in (30, 64, 100, 200):
            b = mp.mpf(beta)
            sig = 2/LAM*mp.quad(lambda v: r(v)*Pcon(v, b), pts)
            absb = 2/LAM*rabs*mp.legenp(-0.5, 0, mp.cosh(b), type=3).real
            F = fcert(beta)
            pr(f'  beta={beta:4d}: signed G_r^low = {mp.nstr(sig, 6):>12}  signed/F_cert = {mp.nstr(sig/F, 4):>10}'
               f'   |abs bound|/F_cert = {mp.nstr(absb/F, 4):>9}   F_cert = {mp.nstr(F, 5)}')
    open(os.path.join(HERE, 'cache', 'tc_lownu_dv.txt'), 'w').write('\n'.join(out) + '\n')
