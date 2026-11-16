"""Referee: independent checks of ADS4_FUNCTIONAL.md (mpmath, no import from adsfunc/).
P_{-1/2+i nu}(cosh beta) = 2F1(1/2-i nu, 1/2+i nu; 1; -sinh^2(beta/2)) (exact hypergeometric), p-integrals by mp.quad.
"""
import sys
import mpmath as mp
from fractions import Fraction as Fr
mp.mp.dps = 20

psi2 = lambda p: (1-p)**3*p*(65+155*p+47*p**2)
psi3 = lambda p: 2*(1-p)**3*(81-1927*p+4601*p**2)

def P(nu, beta):
    z = mp.sinh(beta/2)**2
    return mp.re(mp.hyp2f1(mp.mpf(1)/2-1j*nu, mp.mpf(1)/2+1j*nu, 1, -z))

def W(beta):
    z = mp.sinh(beta/2)**2
    return mp.diff(lambda x: -P(mp.sqrt(x), beta) if x > 0 else -P(0, beta), mp.mpf('1e-12'))  # crude; replaced below

def W_series(beta):
    # -d/d(nu^2) of sum_n (1/2-i nu)_n (1/2+i nu)_n /(n!)^2 (-z)^n at nu=0
    z = mp.sinh(beta/2)**2
    # (1/2-i nu)_n(1/2+i nu)_n = prod_{j<n} ((j+1/2)^2 + nu^2); derivative at nu^2=0: prod * sum 1/(j+1/2)^2
    s = mp.mpf(0); n = 0
    while True:
        prod = mp.mpf(1); acc = mp.mpf(0)
        for j in range(n):
            prod *= (j+mp.mpf(1)/2)**2; acc += 1/(j+mp.mpf(1)/2)**2
        t = prod*acc/mp.factorial(n)**2*(-z)**n
        s += t; n += 1
        if n > 5 and abs(t) < mp.mpf(10)**(-18) * (abs(s)+1): break
        if n > 400: break
    return -s

def FR_parts(h, Lam, m, beta):
    sig = 1 if h == '++' else -1
    u = (Lam/m)**2
    I2 = mp.quad(lambda p: psi2(p)*P(Lam*p, beta), [0, 0.25, 0.5, 0.75, 1])
    I3 = mp.quad(lambda p: psi3(p)*P(Lam*p, beta), [0, 0.25, 0.5, 0.75, 1])
    B2 = 2*I2; B3 = sig*u*I3
    B4 = 10*W_series(beta)*Lam**2*u**2
    B5 = 5*sig*u**3*P(0, beta)
    return B2, B3, B4, B5

def flat_parts(h, J, mfl):
    sig = 1 if h == '++' else -1
    u = 1/mfl**2
    if h == '++':
        D = lambda x: mp.legendre(J, x)
        Phi = 10*J*(J+1)
    else:
        D = lambda x: mp.jacobi(J-4, 0, 8, x)
        Phi = 10*J*(J+1) - 210
    B2 = mp.quad(lambda p: psi2(p)*(2-u*p**2)*D(1-2*u*p**2), mp.linspace(0, 1, 9))
    B3 = sig*u*mp.quad(lambda p: psi3(p)*D(1-2*u*p**2), mp.linspace(0, 1, 9))
    return B2, B3, Phi*u**3

def kin(tau, J):
    al = tau-2; ab = tau+2*J-1
    return mp.sqrt(al*ab), mp.acosh((al**2+ab**2)/(2*al*ab))

if __name__ == '__main__':
    what = sys.argv[1]
    if what == 'exact':
        # J=4, +-: d~^4_44 = 1 -> FF(4,u) = 2 int psi2 - u int p^2 psi2 - u int psi3 - 10 u^3
        from sympy import symbols, integrate, Rational, nsimplify
        p, u = symbols('p u')
        s2 = (1-p)**3*p*(65+155*p+47*p**2); s3 = 2*(1-p)**3*(81-1927*p+4601*p**2)
        a = 2*integrate(s2, (p, 0, 1)); b = integrate(p**2*s2, (p, 0, 1)) + integrate(s3, (p, 0, 1))
        print('2 int psi2 =', a, '  int p^2psi2 + int psi3 =', b)
        for Dg in [100, 1000]:
            uu = Rational(Dg-1, Dg+7)
            v = a - b*uu - 10*uu**3
            print('Dg=%d u_max=%s FF=%.5f  linear approx 0.0603+258.2/(Dg+7)=%.5f' % (Dg, uu, float(v), 19/315+(41/18+30)*8/(Dg+7)))
        # L_AdS asymptotics
        for Lam in [100, 1000]:
            L = mp.quad(lambda q: psi2(q)/(q**2+mp.mpf(9)/(4*Lam**2)), [0, 1/mp.mpf(Lam), 0.1, 1])
            print('Lam=%d L_AdS=%.6f  65log(2Lam/3)-1171/15=%.6f' % (Lam, L, 65*mp.log(2*mp.mpf(Lam)/3)-mp.mpf(1171)/15))
        # W expansion check
        for bt in [0.1, 0.5]:
            z = mp.sinh(bt/2)**2
            print('beta=%.2f W=%.10f  z-5z^2/8=%.10f' % (bt, W_series(bt), z-5*z**2/8))
    elif what == 'island':
        Dg = int(sys.argv[2]); Lam = mp.sqrt((Dg-2)*(Dg-1))
        kap = mp.mpf(10)**0.8
        best = None
        for Lb in [7.5, 8.0, 8.5, 8.7, 9.0, 9.5]:
            beta = Lb/Lam
            pr = FR_parts('++', Lam, kap*Lam, beta)
            F = sum(pr); rel = F/sum(abs(x) for x in pr)
            # flat at same (u, bhat): J from bhat=(2J+1)/mfl
            bh = 2*Lam*mp.sinh(beta/2)
            print('Dg=%d m/Lam=%.4f Lam*beta=%.2f  FR=%+.4e rel=%+.4e  parts=%s  bhat=%.3f' % (Dg, kap, Lb, F, rel, [float(x) for x in pr], bh))
        # flat ++ on nearby physical J (even) at m_fl = kap
        for J in [24, 26, 28, 30]:
            fp = flat_parts('++', J, kap)
            print('  flat ++ J=%d m_fl=%.3f bhat=%.3f m^4F=%+.4e parts=%s' % (J, kap, (2*J+1)/kap, sum(fp), [float(x) for x in fp]))
    elif what == 'nearzero':
        # curvature shift (hybrid) at m_fl=21.7, J=244, +-, for given Dg; both forward-term variants
        Dg = int(sys.argv[2]); Lam = mp.sqrt((Dg-2)*(Dg-1)); mfl = mp.mpf('21.7'); J = 244
        m = mfl*Lam; u = 1/mfl**2; bh = (2*J+1)/mfl
        # tau from m, J
        Bq = 2*J-3; Cq = -2*(2*J-1)-m*m; tau = (-Bq+mp.sqrt(Bq*Bq-4*Cq))/2
        _, beta = kin(tau, J)
        dK = lambda p: P(Lam*p, beta) - mp.besselj(0, bh*p)
        pts = mp.linspace(0, 1, 41)
        dB2 = mp.quad(lambda p: psi2(p)*(2-u*p**2)*dK(p), pts)
        dB3 = -u*mp.quad(lambda p: psi3(p)*dK(p), pts)
        P0 = P(0, beta)
        fw_common = u**3*(-mp.mpf(5)/2 + 10*m*m*W_series(beta) - 10*J*(J+1))
        fw_doc = 5*2*(1-P0)*u**3            # 5(1-sigma)(1-P0), sigma=-1  (as in ADS4_FUNCTIONAL.md §2.4 / common.hybrid_shift)
        fw_ref = 5*(-1)*(P0-1)*u**3          # Regge minus its flat limit for the B5 term: 5 sigma (P0 - 1)
        tot_doc = dB2+dB3+fw_common+fw_doc; tot_ref = dB2+dB3+fw_common+fw_ref
        fl = flat_parts('+-', J, mfl)
        print('Dg=%d Lam=%.2f beta=%.5f  dB2=%.4e dB3=%.4e dFW=%.4e  shift(doc)=%.5e (x Lam^2=%.4f)  shift(ref)=%.5e  B5-term diff=%.2e  flat=%.5e'
              % (Dg, Lam, beta, dB2, dB3, fw_common+fw_doc, tot_doc, tot_doc*Lam**2, tot_ref, fw_doc-fw_ref, sum(fl)))
    elif what == 'j4':
        # J=4 threshold +- curvature shift, doc vs corrected B5 term
        for Dg in [20, 100, 1000]:
            Lam = mp.sqrt((Dg-2)*(Dg-1)); m, beta = kin(mp.mpf(Dg), 4); u = (Lam/m)**2
            bh = 9/(m/Lam)
            dK = lambda p: P(Lam*p, beta) - mp.besselj(0, bh*p)
            pts = mp.linspace(0, 1, 9)
            dB2 = mp.quad(lambda p: psi2(p)*(2-u*p**2)*dK(p), pts)
            dB3 = -u*mp.quad(lambda p: psi3(p)*dK(p), pts)
            P0 = P(0, beta)
            fwc = u**3*(-mp.mpf(5)/2 + 10*m*m*W_series(beta) - 200)
            sd = dB2+dB3+fwc+10*(1-P0)*u**3; sr = dB2+dB3+fwc+5*(1-P0)*u**3
            print('Dg=%d u=%.5f shift(doc)=%+.5e shift(ref)=%+.5e  (x Lam^2: %.1f, %.1f)' % (Dg, u, sd, sr, sd*Lam**2, sr*Lam**2))
