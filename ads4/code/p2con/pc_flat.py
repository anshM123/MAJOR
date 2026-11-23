"""p2con/pc_flat.py (P2_CONTRACT.md §2.1): flat scaling limit of CHMRS's scalar Regge-moment series at O(1/m^2) and O(1/m^4)
(ancillary regge_moments_order8.m, read via p2spin/s5_order4.chmrs_order2, Legendre-reduced as p2spin/s5b_dphi) versus the y-expansion of
the flat partial-wave image E_2(y) d^J(1-2y) (E_2 = (2-y)/(1-y)^2) and of the hybrid floor's stripped (2-y) Dtilde, at fixed b.
Scaling: P -> J0(b), P'(eta) -> -(nu^2/b) J1(b), eta^2 - 1 -> b^2/nu^2, nu^2 = y m^2.  -> cache/flat.txt
Usage: python3 pc_flat.py <path to regge_moments_order8.m>"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'p2spin'))
import numpy as np, sympy as sp
from scipy.special import j0, j1
import s5_order4 as S5, s3_mech as M
import pc_p4 as P4
out = []
pr = lambda s: (print(s, flush=True), out.append(s))
w, bsym = sp.symbols('w b', positive=True)
J0s, J1s = sp.symbols('J0 J1')


def reduce(t2, et):
    sol = {0: (1, 0), 1: (0, 1)}
    for k in range(0, 8):
        c0 = -(2*(k + 1)*et*sol[k + 1][0] + (k*(k + 1) + M.nu2 + sp.Rational(1, 4))*sol[k][0])/(et**2 - 1)
        c1 = -(2*(k + 1)*et*sol[k + 1][1] + (k*(k + 1) + M.nu2 + sp.Rational(1, 4))*sol[k][1])/(et**2 - 1)
        sol[k + 2] = (c0, c1)
    r0 = sp.cancel(sum(t2.coeff(sp.Symbol(f'G{k}'))*sol[k][0] for k in range(9)))
    r1 = sp.cancel(sum(t2.coeff(sp.Symbol(f'G{k}'))*sol[k][1] for k in range(9)))
    return r0, r1


def flat_part(r0, r1, et, order):
    """sum of the flat-scaling terms at O(y^order): returns sympy expression in J0, J1, b."""
    tot = 0
    for r, kind in ((r0, 'P'), (r1, 'D')):
        e = sp.expand(sp.series(sp.cancel(r.subs(et, sp.sqrt(1 + w))), w, 0, 6).removeO())
        for term in sp.Add.make_args(e):
            pw = term.as_powers_dict(); i = int(pw.get(w, 0)); jn = int(pw.get(M.nu2, 0))
            c = term/(w**i*M.nu2**jn)
            if kind == 'P' and 2*jn - 2*i == 2*order:
                tot += c*bsym**(2*i)*J0s
            if kind == 'D' and 2*jn - 2*i + 2 == 2*order:
                tot += -c*bsym**(2*i - 1)*J1s
    return sp.expand(tot)


def main():
    path = sys.argv[1]
    t2, nu, et = S5.chmrs_order2(path, Dphi=3, alpha=7)
    t2 = sp.expand(t2.subs(nu, sp.sqrt(M.nu2)))
    r0, r1 = reduce(t2, et)
    g2 = flat_part(r0, r1, et, 2)
    pr(f'CHMRS O(1/m^4) (Delta_phi = 3) flat-scaling part at O(y^2):  g2(b) = {g2}')
    # O(1/m^2) from the closed form (s3_mech.chmrs_seriesPi, alpha = 7), for the record
    cP, cD = M.chmrs_seriesPi(alpha=7)
    g1 = flat_part(sp.expand(2*cP.subs(M.eta, et)), sp.expand(2*cD.subs(M.eta, et)), et, 1) if False else None
    g2f = sp.lambdify((bsym, J0s, J1s), g2)
    for name, env in (('E_2 d^J (AdS k=2 flat limit)', lambda y: (1 - y)**-2), ('stripped (2-y) D~ (hybrid floor)', lambda y: 1.0)):
        bs = np.linspace(0.6, 4.2, 10)
        dev = [P4.expand_flat('++', b, env)[1] - g2f(b, j0(b), j1(b)) for b in bs]
        pr(f'   ++, {name}: max |c2_flat - g2| over b in [0.6, 4.2] = {max(abs(d) for d in dev):.2e}; at b = 2.7: c2 = '
           f'{P4.expand_flat("++", 2.7, env)[1]:+.6f}, g2 = {g2f(2.7, j0(2.7), j1(2.7)):+.6f}')
    c1 = lambda kc, sh: (lambda b: sh*j0(b) + (11/6)*j0(b) - b*j1(b)/6 - kc*j1(b)/b)
    pr('O(y^2) with mpmath extraction (c1 removed exactly: ++ kC = 1/6, +- kC = -31.8339 (pc_run, beta -> 0); hybrid ++: -2 J0, +-: +2 J0):')
    check2(g2f, {('++', 'E_2 d^J'): c1(1/6, 0), ('++', 'hybrid stripped'): c1(1/6, -2), ('+-', 'E_2 d^J'): c1(-31.8339, 0),
                 ('+-', 'hybrid stripped'): c1(-31.8339, 2)})
    open(os.path.join(HERE, 'cache', 'flat.txt'), 'w').write('\n'.join(out) + '\n')


def c2_extract(h, b, env, c1):
    """O(y^2) coefficient of the flat image at fixed b, after removing J0 + c1 y (mpmath, 40 digits; J = 12..260)."""
    import mpmath as mp
    mp.mp.dps = 40
    Js = list(range(12, 261, 8))
    ys, vals = [], []
    for J in Js:
        y = (mp.mpf(b)/(2*J + 1))**2
        x = 1 - 2*y
        D = mp.legendre(J, x) if h == '++' else mp.jacobi(J - 4, 0, 8, x)*(1 - y)**4
        v = (2 - y)*env(y)*D/2 - mp.besselj(0, b) - c1*y
        ys.append(y); vals.append(v/y**2)
    y0 = max(ys)
    A = mp.matrix([[(ys[i]/y0)**k for k in range(7)] for i in range(len(ys))]); bb = mp.matrix(vals)
    sol = mp.qr_solve(A, bb)[0]
    return float(sol[0])


def check2(g2f, c1map):
    import mpmath as mp
    for h in ('++', '+-'):
        for name, env in (('E_2 d^J', lambda y: (1 - y)**-2), ('hybrid stripped', (lambda y: 1) if h == '++' else (lambda y: (1 - y)**-4))):
            devs = []
            for b in (0.8, 1.7, 2.7, 3.6):
                c2 = c2_extract(h, b, env, c1map[(h, name)](b))
                devs.append((b, c2, g2f(b, j0(b), j1(b))))
            pr(f'   {h} {name:16s}: (b, c2_flat, g2_AdS-scalar) = ' + '; '.join(f'({b}, {c:+.6f}, {g:+.6f})' for b, c, g in devs))


if __name__ == '__main__':
    main()
