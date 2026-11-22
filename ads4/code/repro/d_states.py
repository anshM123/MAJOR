"""repro/d_states.py -- Result D: per-state values of Omega^EE in the hybrid near-gap model (flat NF02 + leading-Regge AdS curvature),
the model that produces the binding twist-edge values (DG_SCAN §2.2 'baseline lambda_0'; EE_DG1000_LARGEJ 'hyb').

Own implementation of (definition read from spinfunc/ee1000_common.hyb and omegau/nf_common.flat_parts_nf; nothing imported):
  B2 = int psi2 (2 - y) Dt,  B3 = int s u psi3 Dt,   y = u p^2, u = (Lam/m)^2,  Dt = P_J(1-2y) (++), P^{(0,8)}_{J-4}(1-2y) (+-);
  cv2 = int psi2 2 (P_{-1/2+i Lam p}(cosh beta) - J0((2J+1) x)),  cv3 = int psi3 u s (P - J0),  x = p sqrt(u),
  cosh beta = 1 + (2J+1)^2/(2 m^2);  lambda/N_s = (B2+B3+cv2+cv3)/(|B2+cv2| + |B3+cv3|)  (diagonal; leading-Regge off-diagonal = 0).
psi2, psi3: NF02 normalised to k(0) = 1, evaluated in factored form (1-p)^5 p X(p)/k0.
P_{-1/2+i nu}(cosh beta): own Mehler-Dirichlet quadrature (sqrt2/pi) int_0^beta cos(nu t) dt/sqrt(cosh beta - cosh t), t = beta sin th;
cross-checked against mpmath.legenp.
NOT reproduced: the beyond-leading-Regge crystal corrections (the -0.0005 that turns 0.2043 into the adopted 0.2038 at 300, and every
corrected value of §5.2/§5.3), which need the near-gap chains.
Usage: python3 d_states.py -> stdout, cache/d_states.txt
"""
import json, os
from fractions import Fraction as Fr
import numpy as np
from scipy.special import j0, roots_legendre
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')
out = []


def pr(s=''):
    print(s, flush=True); out.append(s)


nf = json.load(open(os.path.join(ROOT, 'robust', 'functionals_noforward.json')))['NF02']
s2 = [Fr(x) for x in nf['psi2_exact']]; s3 = [Fr(x) for x in nf['psi3_exact']]
mom = lambda c, j: sum(cn/(n + j + 1) for n, cn in enumerate(c))
k0 = float(-2*mom(s2, 6) + 2*mom(s3, 4))
X = np.array([float(Fr(x)) for x in nf['X']])/k0; Y = np.array([float(Fr(x)) for x in nf['Y']])/k0
psi2 = lambda p: (1 - p)**5*p*np.polynomial.polynomial.polyval(p, X)
psi3 = lambda p: (1 - p)**5*p*np.polynomial.polynomial.polyval(p, Y)


def gl01(n):
    x, w = roots_legendre(n); return (x + 1)/2, w/2


TH, WTH = roots_legendre(400); TH = (TH + 1)*np.pi/4; WTH = WTH*np.pi/4


def legP_mehler(nu, beta):
    """P_{-1/2+i nu}(cosh beta), real nu array."""
    t = beta*np.sin(TH)
    den = np.sqrt(np.maximum(np.cosh(beta) - np.cosh(t), 0))
    # stable: cosh b - cosh t = 2 sinh((b+t)/2) sinh((b-t)/2)
    den = np.sqrt(2*np.sinh((beta + t)/2)*np.sinh((beta - t)/2))
    ker = beta*np.cos(TH)/den*WTH
    return (np.sqrt(2)/np.pi)*np.cos(np.outer(nu, t)) @ ker


def jac_one(n, al, be, x):
    p0 = np.ones_like(x)
    if n == 0:
        return p0
    p1 = (al - be)/2 + (al + be + 2)*x/2
    for k in range(2, n + 1):
        A = 2*k*(k + al + be)*(2*k + al + be - 2)
        B = (2*k + al + be - 1)*((2*k + al + be)*(2*k + al + be - 2)*x + al*al - be*be)
        C = 2*(k + al - 1)*(k + be - 1)*(2*k + al + be)
        p0, p1 = p1, (B*p1 - C*p0)/A
    return p1


def state(Dg, J, mfl, h, npn=600):
    Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0)); u = 1/mfl**2; s = 1 if h == '++' else -1
    p, w = gl01(npn); y = u*p*p; x = p*np.sqrt(u)
    Dt = jac_one(J, 0, 0, 1 - 2*y) if h == '++' else jac_one(J - 4, 0, 8, 1 - 2*y)
    B2 = np.sum(w*psi2(p)*(2 - y)*Dt); B3 = np.sum(w*s*u*psi3(p)*Dt)
    m = mfl*Lam; beta = np.arccosh(1 + (2*J + 1)**2/(2*m*m))
    P = legP_mehler(Lam*p, beta); jz = j0((2*J + 1)*x)
    cv2 = np.sum(w*psi2(p)*2*(P - jz)); cv3 = np.sum(w*psi3(p)*u*s*(P - jz))
    Ns = abs(B2 + cv2) + abs(B3 + cv3)
    return dict(val=(B2 + B3 + cv2 + cv3)/Ns, flat=(B2 + B3)/(abs(B2) + abs(B3)), B2=B2, B3=B3, cv2=cv2, cv3=cv3)


def mgap(Dg, J):
    al = Dg - 2.0; return np.sqrt(al*(al + 2*J + 1))/np.sqrt((Dg - 2.0)*(Dg - 1.0))


if __name__ == '__main__':
    pr('# repro/d_states.py: hybrid (flat NF02 + leading-Regge curvature) per-state values lambda/N_s')
    # Legendre check
    for nu, eta in ((50.0, 1.0007), (250.0, 1.00066), (3.0, 1.05)):
        mine = legP_mehler(np.array([nu]), np.arccosh(eta))[0]
        ref = float(mp.re(mp.legenp(-0.5 + 1j*nu, 0, eta, type=3)))
        pr(f'  Legendre check P_(-1/2+{nu}i)({eta}): Mehler {mine:+.12e}  mpmath {ref:+.12e}')
    for Dg in (300, 1000):
        pr(f'\nDg = {Dg}: twist-edge masses mgap(J=4,5) = {mgap(Dg, 4):.5f}, {mgap(Dg, 5):.5f} Lam')
        J = 5
        for mfl in sorted({mgap(Dg, 5), 1.0048, 1.01, 1.02, 1.05}):
            if mfl < mgap(Dg, 5) - 1e-9:
                tag = '  (below the twist edge: not a physical state)'
            else:
                tag = ''
            st = state(Dg, J, mfl, '+-')
            pr(f'  J = 5 +-, m = {mfl:.5f} Lam: hybrid {st["val"]:.5f}  (flat NF02 alone {st["flat"]:.5f}; cv2 {st["cv2"]:+.3e}, cv3 {st["cv3"]:+.3e}){tag}')
        # per-J minimum over masses twist edge .. 1.3 Lam (hybrid only), for comparison with the CL/HC adopted (corrected) values
        for J in range(0, 17):
            vals = []
            for h in ('++', '+-'):
                if (h == '++' and J % 2) or (h == '+-' and J < 4):
                    continue
                for mfl in np.concatenate([[mgap(Dg, J)], np.linspace(max(1.0, mgap(Dg, J)) + 0.005, 1.3, 30), [1.5, 2, 3, 4, 5]]):
                    vals.append((state(Dg, J, mfl, h)['val'], mfl, h))
            if not vals:
                continue
            v = min(vals)
            pr(f'  J = {J:2d}: hybrid min over m in [edge, 5 Lam] = {v[0]:+.4f} at {v[1]:.4f} Lam ({v[2]})')
    # flat NF02 margin on a grid (certificate claims min F/N = 0.200)
    pr('\nflat NF02: min (B2+B3)/(|B2|+|B3|) over a grid of (u, J):')
    for h in ('++', '+-'):
        best = (9,)
        for u in np.unique(np.concatenate([np.logspace(-3.5, 0, 90), np.linspace(0.85, 1, 31)])):
            Jmax = int(min(300, max(12, 40/np.sqrt(u))))
            p, w = gl01(Jmax + 40); y = u*p*p; s = 1 if h == '++' else -1
            al, be, J0 = (0, 0, 0) if h == '++' else (0, 8, 4)
            for J in range(J0, Jmax + 1, 1 if h == '+-' else 2):
                Dt = jac_one(J - J0, al, be, 1 - 2*y)
                B2 = np.sum(w*psi2(p)*(2 - y)*Dt); B3 = np.sum(w*s*u*psi3(p)*Dt)
                r = (B2 + B3)/(abs(B2) + abs(B3))
                if r < best[0]:
                    best = (r, u, J)
        pr(f'  {h}: min {best[0]:.5f} at u = {best[1]:.4f} (m = {1/np.sqrt(best[1]):.3f}), J = {best[2]}')
    open(os.path.join(HERE, 'cache', 'd_states.txt'), 'w').write('\n'.join(out) + '\n')
