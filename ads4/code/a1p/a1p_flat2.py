"""
a1p_flat2 (A1PRIME.md section 4.3-4.4): exact sign certificates (sympy real-root isolation on x = m^2 in [0,1]) for the flat light-matter
actions of a1p_flat, the crossing mass m0 of Q's scalar action, and the repaired upper end Omega^Q_up + eps*Omega^EE (eps = smallest value
making the light-scalar action non-negative on all m in [0,1]).  Units M = 1 (flat) <-> Lambda (AdS).
"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
import sympy as sp
import numpy as np
from scipy.integrate import quad
import a1p_flat as A

OUT = os.path.join(HERE, 'cache', 'a1p_flat2.txt')

def main():
    L = []
    pr = lambda s: (print(s, flush=True), L.append(s))
    x = sp.symbols('x')
    d = json.load(open(os.path.join(A.DATA, 'nf02_psi.json')))
    s2 = [Fr(v) for v in d['psi2_exact']]; s3 = [Fr(v) for v in d['psi3_exact']]
    k0 = -2*A.mom(s2, 6) + 2*A.mom(s3, 4)
    muD2 = lambda j: A.mom(s2, j)/k0; muD3 = lambda j: A.mom(s3, j)/k0
    FD0, FD2 = A.matter_polys(muD2, muD3)
    q = json.load(open(os.path.join(A.DATA, 'q3_full_delta0.2.json')))
    q2 = [Fr(float(v)) for v in q['psi2']]; q3 = [Fr(float(v)) for v in q['psi3']]
    a = Fr(float(q['a']))
    muQ2 = lambda j: A.mom(q2, j); muQ3 = lambda j: A.mom(q3, j)
    FQ0, FQ2 = A.matter_polys(muQ2, muQ3, a=a)
    tosp = lambda c: sum(sp.Rational(ci.numerator, ci.denominator)*x**i for i, ci in enumerate(c))
    for nm, c in (('F_D J=0', FD0), ('F_D J=2', FD2), ('F_Q J=0', FQ0), ('F_Q J=2', FQ2)):
        P = sp.Poly(tosp(c), x)
        roots = [float(r) for r in P.real_roots() if 0 <= float(r) <= 1]
        v0 = float(P.eval(0)); v1 = float(P.eval(1))
        pr(f'{nm}: real roots of the exact polynomial with x = m^2 in [0,1]: '
           f'{roots or "none"}; P(0) = {v0:+.6f}, P(1) = {v1:+.6f}'
           f' -> {"sign-definite on [0,1] [P]" if not roots else "changes sign"}')
        if roots:
            r = roots[0]
            pr(f'   crossing at m0^2 = {r:.6f}, m0 = {r**0.5:.6f} M; negative for m < m0')
    # minimal eps with FQ0 + eps FD0 >= 0 on [0,1]: max_x of -FQ0/FD0
    xs = np.linspace(0, 1, 100001)
    ratio = -np.array([A.peval(FQ0, xx) for xx in xs])/np.array([A.peval(FD0, xx) for xx in xs])
    eps = ratio.max(); xa = xs[ratio.argmax()]
    pr(f'\nRepair: eps* = max_m [-F_Q0/F_D0] = {eps:.6f} at m^2 = {xa:.4f} (exact value at m = 0: '
       f'{float(-FQ0[0]/FD0[0]):.8f})')
    eps = float(-FQ0[0]/FD0[0])
    # combined J = 2 action stays positive
    comb2 = [FQ2[i] + Fr(eps)*FD2[i] for i in range(4)]
    P = sp.Poly(tosp(comb2), x)
    pr(f'combined J = 2 action: roots in [0,1]: {[float(r) for r in P.real_roots() if 0 <= float(r) <= 1] or "none"};'
       f' min on grid {min(A.peval(comb2, xx) for xx in xs):+.5f}')
    # AdS constants: L_A = int psi2/(p^2 + 9/(4 Lam^2)) with NF02 (UT: R0 = 240 sqrt(16 L_A/K)/Lam^4, 16 L_A = 1.07e4 / 1.43e4)
    psi2 = lambda p: float(sum(float(c)*p**i for i, c in enumerate(s2))/float(k0))
    psi2q = lambda p: float(sum(float(c)*p**i for i, c in enumerate(q2)))
    ID = quad(psi2, 0, 1)[0]; IQ = quad(psi2q, 0, 1)[0]
    kap = float(q['kappa'])
    pr(f'I_D = int psi2^NF = {ID:.6f} (UT: 10.385947); I_Q = {IQ:.7f} (UT: 1.0160366)')
    for Dg, twoLQ in ((300, 72.041), (1000, 92.814)):
        Lam = np.sqrt((Dg - 2)*(Dg - 1))
        LA = quad(lambda p: psi2(p)/(p*p + 9/(4*Lam*Lam)), 0, 1, limit=400, points=[1e-4, 1e-3, 1e-2])[0]
        R0 = 240*np.sqrt(16*LA)/Lam**4
        up = twoLQ + 2*eps*LA
        cc = 57600*(IQ + eps*ID)**2/((kap + eps/8)*Lam**4)
        pr(f'Delta_g = {Dg}: L_A = {LA:.3f} (16 L_A = {16*LA:.1f}; R0 check {R0:.5e}); published upper end 2L_AdS = {twoLQ}; '
           f'repaired (Omega^Q_up + {eps:.5f} Omega^EE): 2L_AdS + 2 eps L_A = {up:.3f} (x{up/twoLQ:.3f}); '
           f'kappa -> kappa + eps/8 = {kap + eps/8:.5f}; c1^2 coefficient {cc:.4e} (was {57600*IQ**2/(kap*Lam**4):.4e})')
    with open(OUT, 'w') as fh:
        fh.write('\n'.join(L) + '\n')

if __name__ == '__main__':
    main()
