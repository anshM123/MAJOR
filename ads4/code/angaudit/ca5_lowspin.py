"""
Audit of ANGULAR_ROUTE §21 item 3, low-spin exceptional sectors (J = 0, 2; J = 1, 3 have no structures).
Collinear F per unit coefficient (from ca1_collinear, O polarised e_+, repo canonical J=0,2 normalisations):
   J=0 even: F_{-2,-2} = (s-4)(s-3)/(16 s(s+1));    J=0 odd: F_{-2,-1} = i(2s-7)/8          (s = D/2)
   J=2 even: F_{-2,-2} = -(s-6)(s-5)(s-4)/8;         J=2 odd: F_{-2,-1} = i(2s-11)(2s-9)/8  (s = D/2 + 1)
Crossing budget per dyadic Delta-shell: sum 4^D 2^-J |F|^2 <~ D^{11} (S ~ U ~ D/2, all M).
We print, per dyadic shell [D0, 2 D0] (unit-spaced operators, budget spread evenly),
   sum Delta^2 |a|^2 = sum Delta^2 |K|^2 |lambda|^2,  |lambda|^2 = budget_per_op / (4^D 2^-J |F/lambda|^2),
and the large-Delta exponents of |K|^2/4^D.
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'kernels'))
import mpmath as mp
import closed_forms as cf
mp.mp.dps = 30

def K2J2e(D):
    return (mp.mpf(9)/2*mp.pi**3.5*mp.mpf(4)**D*(D - 3)*(D - 2)*mp.gamma(D - 0.5)
            *mp.rgamma(D + 2)*mp.rgamma(2 + D/2)**2*mp.rgamma(3 - D/2)**2)
def K2J2o(D):
    return (mp.mpf(9)/8*mp.pi**4.5*(D - 3)*(D - 2)*mp.gamma(D + 2)*mp.gamma(D - 0.5)
            *mp.rgamma(2 + D/2)**2*mp.rgamma((D + 3)/2)**4*mp.rgamma((7 - D)/2)**2)
# cross-check the J=2 closed forms against cf amplitude functions
for D in (mp.mpf('11.3'), mp.mpf('57.9')):
    e = cf.K2(cf.M0_J2_even(D), D, 2, 0)/K2J2e(D) - 1; o = cf.K2(cf.M0_J2_odd(D), D, 2, 0)/K2J2o(D) - 1
    print(f"J=2 closed-form K^2 consistency at D={D}: {mp.nstr(e, 3)}, {mp.nstr(o, 3)}")

sectors = {
    'J0 even': (0, cf.K0_even_scalar2, lambda D: (D/2 - 4)*(D/2 - 3)/(16*(D/2)*(D/2 + 1))),
    'J0 odd ': (0, cf.K0_odd_scalar2, lambda D: (D - 7)/8),
    'J2 even': (2, K2J2e, lambda D: (D/2 + 1 - 6)*(D/2 + 1 - 5)*(D/2 + 1 - 4)/8),
    'J2 odd ': (2, K2J2o, lambda D: (D + 2 - 11)*(D + 2 - 9)/8),
}
for name, (J, K2, Ff) in sectors.items():
    # exponent of |K|^2/4^D
    D1, D2 = mp.mpf(2000) + mp.mpf(0.5), mp.mpf(4000) + mp.mpf(0.5)
    ex = mp.log((K2(D2)/mp.mpf(4)**D2)/(K2(D1)/mp.mpf(4)**D1))/mp.log(D2/D1)
    out = []
    for D0 in (32, 128, 512, 2048):
        tot = mp.mpf(0)
        for j in range(D0):
            D = mp.mpf(D0 + j) + mp.mpf('0.5')
            lam2 = (mp.mpf(D0)**11/D0)/(mp.mpf(4)**D*mp.mpf(2)**(-J)*abs(Ff(D))**2)
            tot += D**2*K2(D)*lam2
        out.append((D0, tot))
    print(f"{name}: |K|^2/4^D ~ D^{mp.nstr(ex, 4)};  per-shell sum D^2|a|^2: " +
          ", ".join(f"D0={d}: {mp.nstr(t, 4)}" for d, t in out) +
          f";  ratio 2048/512 = {mp.nstr(out[3][1]/out[2][1], 4)} (D^-3/2 -> {mp.nstr(mp.mpf(4)**-1.5, 4)})")
