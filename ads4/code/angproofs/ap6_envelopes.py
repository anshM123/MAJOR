"""Step 6: factorised absolute kernels and envelope bounds with explicit constants (Wendel's inequality only).
(i)  high-precision check that the factorised forms equal the KERNELS_D3 closed forms (closed_forms.py);
(ii) explicit constants: |K_X|^2 <= c_X(J, Delta_0) * 4^Delta 2^{-J} sqrt(J+1) (sh)^{+-1/2} Delta^{-15/2}, all J >= 4,
     h >= 5 (Delta >= 14);  (iii) numerical scan of the ratio."""
import mpmath as mp
from common import *
import closed_forms as cf
mp.mp.dps = 40

def brackets(Delta, J):
    s = (Delta + J)/2; h = (Delta - J)/2
    cJ = mp.mpf(2)**J/mp.binomial(2*J, J)
    gs = mp.gamma(s + 0.5)/mp.gamma(s); gh = mp.gamma(h)/mp.gamma(h - 0.5)
    return s, h, cJ, gs, gh

def K2_fact(name, Delta, J):
    s, h, cJ, gs, gh = brackets(Delta, J)
    D0 = mp.gamma(Delta - 1)*mp.gamma(Delta - 0.5)/mp.gamma(Delta + 3)**2
    D4 = mp.gamma(Delta - 5)*mp.gamma(Delta - 0.5)/(mp.gamma(Delta + 3)*mp.gamma(Delta - 1))
    rho = mp.binomial(2*J, J)/mp.binomial(2*J, J + 4) if J >= 4 else 0
    sn2, cs2 = mp.sin(mp.pi*h)**2, mp.cos(mp.pi*h)**2
    p32 = mp.pi**1.5*4**mp.mpf(Delta)
    return {'pp': 864*p32*sn2*cJ*gs*gh*D0,
            'odd': 3456*p32*cs2*cJ/gs/gh*D0,
            'pm': 216*p32*sn2*cJ*rho*gs*gh*D4,
            'oddJ': 3456*p32*cs2*cJ*rho/gs/gh*D4}[name]

def K2_closed(name, Delta, J):
    return {'pp': cf.Kpp2, 'odd': cf.Kodd2, 'pm': cf.Kpm2, 'oddJ': cf.KoddJ2}[name](Delta, J)

def envelope(name, Delta, J):
    s = (Delta + J)/mp.mpf(2); h = (Delta - J)/mp.mpf(2)
    e = 4**mp.mpf(Delta)*mp.mpf(2)**(-J)*mp.sqrt(J + 1)/mp.mpf(Delta)**7.5
    return e*mp.sqrt(s*h) if name in ('pp', 'pm') else e/mp.sqrt(s*h)

# explicit constants.  Wendel (1948): for x > 0, 0 < a < 1:  (x/(x+a))^{1-a} <= Gamma(x+a)/(x^a Gamma(x)) <= 1.
#   (W1) Gamma(x+1/2)/Gamma(x) <= sqrt(x);   (W2) Gamma(x+1/2)/Gamma(x) >= x/sqrt(x+1/2).
#   2^J/C(2J,J) = 2^{-J} sqrt(pi) Gamma(J+1)/Gamma(J+1/2) <= 2^{-J} sqrt(pi) sqrt(J+1/2)            (W1, x = J+1/2)
#   Gamma(s+1/2)/Gamma(s) <= sqrt(s);  Gamma(h)/Gamma(h-1/2) <= sqrt(h-1/2) <= sqrt(h)               (W1)
#   Gamma(s)/Gamma(s+1/2) <= sqrt(s+1/2)/s;  Gamma(h-1/2)/Gamma(h) <= sqrt(h)/(h-1/2)                 (W2)
#   Gamma(Delta+5/2)/Gamma(Delta+3) <= sqrt(Delta+3)/(Delta+5/2)                                      (W2)
def g0(D):
    """sup_{Delta >= D} Delta^{15/2} Gamma(Delta-1)Gamma(Delta-1/2)/Gamma(Delta+3)^2.  Using (W2),
    Gamma(Delta-1/2)/Gamma(Delta+3) <= sqrt(Delta+3)/((Delta+5/2)(Delta-1/2)(Delta+1/2)(Delta+3/2)) and
    sqrt(Delta(Delta+3)) <= Delta+5/2, every factor Delta/(Delta+k) (k >= 0) is <= 1, leaving
    Delta/(Delta-1) * Delta/(Delta-1/2), decreasing in Delta."""
    return D/(D - 1)*D/(D - mp.mpf(1)/2)
def g4(D):
    """sup_{Delta >= D} Delta^{15/2} Gamma(Delta-5)Gamma(Delta-1/2)/(Gamma(Delta+3)Gamma(Delta-1)), same method"""
    return D**4/((D - 5)*(D - 4)*(D - 3)*(D - 2))*D/(D - mp.mpf(1)/2)

def constants(D0=14, s0=9, h0=5, J0=4):
    sq = mp.sqrt(mp.pi)            # from 2^J/C(2J,J) <= 2^{-J} sqrt(pi) sqrt(J+1)
    rho0 = mp.binomial(2*J0, J0)/mp.binomial(2*J0, J0 + 4)
    odd_fac = mp.sqrt(1 + 1/(2*mp.mpf(s0)))*h0/(h0 - 0.5)
    return {'pp': 864*mp.pi**1.5*sq*g0(mp.mpf(D0)),
            'odd': 3456*mp.pi**1.5*sq*odd_fac*g0(mp.mpf(D0)),
            'pm': 216*mp.pi**1.5*sq*rho0*g4(mp.mpf(D0)),
            'oddJ': 3456*mp.pi**1.5*sq*rho0*odd_fac*g4(mp.mpf(D0))}

if __name__ == '__main__':
    print("(i) factorised form vs KERNELS_D3 closed forms (max rel. deviation):")
    for name, Js in (('pp', [4, 6, 10, 30, 200]), ('odd', [4, 8, 30, 200]), ('pm', [4, 6, 30, 200]), ('oddJ', [5, 9, 31, 201])):
        dev = 0
        for Jv in Js:
            for h in (mp.mpf('5.25'), mp.mpf('7.6'), mp.mpf('40.3'), mp.mpf('300.1')):
                Dv = Jv + 2*h
                dev = max(dev, abs(K2_fact(name, Dv, Jv)/K2_closed(name, Dv, Jv) - 1))
        print(f"   {name:5s}: {mp.nstr(dev, 3)}")
    C = constants()
    print("(ii) explicit constants c_X (|K_X|^2 <= c_X * envelope for all J >= 4, h >= 5):")
    for k, v in C.items():
        print(f"   {k:5s}: {mp.nstr(v, 6)}")
    print("(iii) scan of |K|^2/envelope (max over grid) vs constant:")
    for name in ('pp', 'odd', 'pm', 'oddJ'):
        mx = 0
        Jlist = [4, 5, 6, 8, 12, 20, 50, 100, 400, 2000, 10000] if name in ('pp', 'pm', 'odd') else [5, 7, 9, 13, 21, 51, 101, 401, 2001, 10001]
        for Jv in Jlist:
            if name in ('pp', 'pm', 'odd') and Jv % 2: continue
            for h in [5 + mp.mpf(k)/8 for k in range(8)] + [mp.mpf(x) + 0.37 for x in (7, 11, 30, 100, 1000, 10000)]:
                Dv = Jv + 2*h
                r = K2_fact(name, Dv, Jv)/envelope(name, Dv, Jv)
                mx = max(mx, r)
        print(f"   {name:5s}: max ratio {mp.nstr(mx, 5)}  <=  c = {mp.nstr(C[name], 5)}: {mx <= C[name]}")
