"""Low-spin sectors J = 0, 2: collinear F (exact in s at fixed J), from the KERNELS_D3 structures (data/*.json),
by the same collinear evaluation as ap3.  O polarised e_+; per unit canonical coefficient."""
import json
from common import *
from ap3_collinear import collinear
DATA = os.path.join(HERE, '..', 'kernels', 'data')
for Jv in (0, 2):
    d = json.load(open(os.path.join(DATA, f'structures_J{Jv}.json')))
    for name, v in d.items():
        vec = {}
        for m, c in zip(v['monomials'], v['coefficients']):
            mr = (m[0], m[1], Jv - m[2], m[3], m[4], m[5], m[6])
            vec[mr] = sp.sympify(c, locals={'s': S})
        F, sgn = collinear(vec, v['parity'])
        F = {k: sp.factor(val.subs(J, Jv)) for k, val in F.items()}
        key = (-2, -2) if v['parity'] == 'even' else (-2, -1)
        Delta = sp.Symbol('Delta')
        print(f"J={Jv} {name:4s} ({v['parity']}): F{key} = {F.get(key, 0)}   [(2V4)^J = {sgn**Jv}]"
              f"   = {sp.factor(F.get(key, 0).subs(S, (Delta + Jv)/2))} in Delta")

# growth check: |K|^2 / |F per unit lambda|^2 ~ 4^Delta Delta^{-q}, q = p - 2g, with (p, g) the kernel and F exponents;
# the per-dyadic-shell budget is Delta^{11}, so sum_shell Delta^2 |a|^2 <~ Delta^{2 + 11 - q - ...}
import mpmath as mp
import closed_forms as cf
mp.mp.dps = 30
def K2_J2(Delta, par):
    D = mp.mpf(Delta)
    if par == 'even':
        return 4.5*mp.pi**3.5*4**D*(D - 3)*(D - 2)*mp.gamma(D - 0.5)/(mp.gamma(D + 2)*mp.gamma(2 + D/2)**2*mp.gamma(3 - D/2)**2)
    return 1.125*mp.pi**4.5*(D - 3)*(D - 2)*mp.gamma(D + 2)*mp.gamma(D - 0.5)/(mp.gamma(2 + D/2)**2*mp.gamma((D + 3)/2)**4*mp.gamma((7 - D)/2)**2)
def K2_J0(Delta, par):
    D = mp.mpf(Delta)
    if par == 'even':
        return 6.75*mp.pi**4.5*mp.gamma(D)*mp.gamma(D - 0.5)/(mp.gamma(2 + D/2)**4*mp.gamma(3 - D/2)**2*mp.gamma((3 + D)/2)**2)
    return 6.75*mp.pi**4.5*mp.gamma(D)*mp.gamma(D - 0.5)/(mp.gamma(2 + D/2)**2*mp.gamma((1 + D)/2)**2*mp.gamma((3 + D)/2)**2*mp.gamma((7 - D)/2)**2)
Fc = {(0, 'even'): lambda D: (D - 8)*(D - 6)/(16*D*(D + 2)), (0, 'odd'): lambda D: (D - 7)/8,
      (2, 'even'): lambda D: (D - 10)*(D - 8)*(D - 6)/64, (2, 'odd'): lambda D: (D - 9)*(D - 7)/8}
print("low spin: sup over the shell [D0, 2D0] of  Delta^{2+11} |K|^2/(4^Delta |F_coef|^2)  (should scale like D0^{-3/2}):")
for (Jv, par), f in Fc.items():
    vals = []
    for D0 in (128, 512, 2048):
        mx = 0
        for q in range(41):
            D = mp.mpf(D0)*(1 + mp.mpf(q)/40) + mp.mpf('0.3')
            K2 = K2_J0(D, par) if Jv == 0 else K2_J2(D, par)
            mx = max(mx, D**13*K2/(4**D*abs(f(D))**2))
        vals.append(mx)
    print(f"   J={Jv} {par:4s}: D0=128,512,2048 -> {[mp.nstr(v, 4) for v in vals]};  ratio per factor 4 in D0: "
          f"{mp.nstr(vals[1]/vals[0], 4)}, {mp.nstr(vals[2]/vals[1], 4)}  (4^-1.5 = 0.125)")
