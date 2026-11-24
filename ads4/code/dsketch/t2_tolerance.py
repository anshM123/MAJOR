"""
t2_tolerance (DSKETCH_UP.md §3.4-3.5): how large an AdS correction to the light-state values the sign statements tolerate.

Exact structure (DSKETCH_UP §3.3): the k = 2, 3 light value of a light J = 0 / J = 2 exchange is a polynomial in nu^2 whose TOP
coefficient is the flat (CHLPS) one and whose lower coefficients carry AdS shifts of O(1) in AdS units (O(Lambda^-2) relative).
In Lambda units (nu = Lambda p, m in Lambda units):
  J = 0:  B2 = -(p^2 - 2 mu^2/Lambda^2),  B3 = 1 (exact)        -> the whole AdS effect is m^2 -> mu^2/Lambda^2, mu^2 in AdS units;
  J = 2:  B2 = -6p^6 + (18 m^2 + s1) p^4 + (-13 m^4 + s2) p^2 + (2 m^6 + s3),  B3 = 6 p^4 + (-6 m^2 + t1) p^2 + (m^4 + t2)
          (per |g2|^2/m^4), shifts bounded as |s_j|, |t_j| <= kap (m^2 + kap)^(d_j - 1), kap = K/Lambda^2, d_j the flat m^2-degree.
Functionals: Omega^EE (NF02 weights, k(0) = 1), Omega^Q_up (Q3-full-delta0.2; J = 2 also gets +12 a from the k = 4 arc), and
the repaired Omega^Q_up + eps Omega^EE with eps = 0.040 and eps* = 0.0380339.
Outputs: J = 0: the critical mu^2 (AdS units) below which each functional turns negative, at Delta_g = 300, 1000;
         J = 2: the largest K (AdS units) for which the worst-case value stays positive for all 0 < m <= Lambda.
Writes cache/t2_tolerance.txt
"""
import os, sys, json
from fractions import Fraction as Fr
sys.dont_write_bytecode = True
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', '..', 'prl', 'data')
lines = []
def log(s=''):
    print(s, flush=True); lines.append(s)

nf = json.load(open(os.path.join(DATA, 'nf02_psi.json')))
k0 = Fr(nf['k0'])
c2D = [Fr(v)/k0 for v in nf['psi2_exact']]; c3D = [Fr(v)/k0 for v in nf['psi3_exact']]
q3 = json.load(open(os.path.join(DATA, 'q3_full_delta0.2.json')))
c2Q = [Fr(v) for v in q3['psi2']]; c3Q = [Fr(v) for v in q3['psi3']]
aQ = Fr(q3['a'])
def mom(c, j):          # int_0^1 psi(p) p^j dp
    return sum(ci/(i + j + 1) for i, ci in enumerate(c))
W = {'D': (c2D, c3D, Fr(0)), 'Q': (c2Q, c3Q, aQ)}
LAM = {300: np.sqrt(298.0*299.0), 1000: np.sqrt(998.0*999.0)}
EPS = 0.040; EPSS = 0.4304302/11.3170270

# ---------------- J = 0
F0 = {}
for t, (c2, c3, a) in W.items():
    F0[t] = (float(mom(c3, 0) - mom(c2, 2)), float(2*mom(c2, 0)))       # F = A + B m^2
log('J = 0 flat polynomials (per |g0|^2, m in Lambda units):')
for t in W:
    log(f'  F_{t}^0 = {F0[t][0]:+.6f} + {F0[t][1]:.6f} m^2')
combos = {'EE': (1.0, 0.0), 'Q_up': (0.0, 1.0), 'Q_up + 0.040 EE': (EPS, 1.0), 'Q_up + eps* EE': (EPSS, 1.0)}
log('J = 0: the AdS light value is F(m^2 -> mu^2/Lambda^2) exactly (B3 exact; B2 constant shifted); critical mu^2 (AdS units):')
for name, (e, q) in combos.items():
    A = e*F0['D'][0] + q*F0['Q'][0]; B = e*F0['D'][1] + q*F0['Q'][1]
    x0 = -A/B
    s = '  ' + f'{name:16s}: F = {A:+.7f} + {B:.6f} mu^2/Lambda^2 ; positive iff mu^2 > {x0:+.6f} Lambda^2 = '
    s += ', '.join(f'{x0*LAM[g]**2:+.1f} (Delta_g = {g})' for g in (300, 1000))
    log(s)

# ---------------- J = 2
def F2_parts(t):
    c2, c3, a = W[t]
    m2 = lambda j: float(mom(c2, j)); m3 = lambda j: float(mom(c3, j))
    # flat polynomial in x = m^2 (per |g2|^2/m^4): F = sum f_i x^i
    f = [-6*m2(6) + 6*m3(4) + 12*float(a), 18*m2(4) - 6*m3(2), -13*m2(2) + m3(0), 2*m2(0)]
    # sensitivities to the shifts: s1 p^4, s2 p^2, s3 (in B2) ; t1 p^2, t2 (in B3)
    sens = {'s1': m2(4), 's2': m2(2), 's3': m2(0), 't1': m3(2), 't2': m3(0)}
    return f, sens
deg = {'s1': 1, 's2': 2, 's3': 3, 't1': 1, 't2': 2}
log('J = 2 flat polynomials (per |g2|^2/m^4):')
for t in W:
    f, _ = F2_parts(t)
    log(f'  F_{t}^2 = ' + ' '.join(f'{v:+.6f} m^{2*i}' for i, v in enumerate(f)))
xs = np.linspace(1e-6, 1.0, 20001)
for name, (e, q) in (('EE', (1.0, 0.0)), ('Q_up', (0.0, 1.0)), ('Q_up + 0.040 EE', (EPS, 1.0))):
    fD, sD = F2_parts('D'); fQ, sQ = F2_parts('Q')
    f = [e*a_ + q*b_ for a_, b_ in zip(fD, fQ)]
    sens = {k: abs(e*sD[k]) + abs(q*sQ[k]) for k in sD}
    flat = sum(v*xs**i for i, v in enumerate(f))
    def worst(kap):
        pen = sum(sens[k]*kap*(xs + kap)**(deg[k] - 1) for k in sens)
        return np.min(flat - pen)
    lo, hi = 0.0, 1.0
    for _ in range(80):
        mid = (lo + hi)/2
        if worst(mid) > 0:
            lo = mid
        else:
            hi = mid
    log(f'  {name:16s}: min flat F^2 = {flat.min():.4f}; positivity survives every shift with K <= {lo:.4e} Lambda^2 = '
        + ', '.join(f'{lo*LAM[g]**2:.3g} (Delta_g = {g})' for g in (300, 1000)) + '  [AdS units]')
# natural size of the shifts: the Barnes map T -> (A^2 + nu^2)/(2A) (t1b_barnes) shifts a Mellin-linear coefficient by A^2
log('Natural size of the AdS shift from the Mellin -> nu map alone (t1b_barnes): nu^2 -> nu^2 + A^2, A = a - 1/2, a = 2*3 + k - 1:')
for k in (2, 3):
    A = 2*3 + k - 1 - 0.5
    log(f'  k = {k}: A = {A}, A^2 = {A*A:.2f}  -> shift of mu^2 by A^2/2 = {A*A/2:.2f} (AdS units)')
open(os.path.join(HERE, 'cache', 't2_tolerance.txt'), 'w').write('\n'.join(lines) + '\n')
