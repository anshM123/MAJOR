"""
v4d_og (RESULT_Q_V4DPRIME.md §2): O(G) u-channel twist-1 exchange at the light-cone-edge/Regge corner, exact light-ray-polarised spinning
structures, raw per-nu and with the actual de-shift weights.

Homogeneity lemma [D]: at the regulated corner (t1 = t3 = 1/2, dphi ~ sqrt(eps)) the invariants r = z1.z3, p = z1.x13, q = z3.x13, X = x13^2
are all O(eps) while x13 itself and every contraction of the exchanged operator with z1, z3, x13 stay O(1).  A leading-twist (13) light-cone
term of a twist-tau exchange between two light-ray-polarised spin-2, Delta = 3 operators is X^{(tau - 6)/2} x (polynomial in r, p, q, X of
z1- and z3-degree 2 each), and dimensional + degree counting gives eps^{-1 + tau/2} for EVERY tensor structure (tau = 0: the identity,
eps^{-1}; tau = 1: T and the twist-1 family, eps^{-1/2}).  Pointwise (eps = 0, dphi -> 0) the same counting gives dphi^{-2 + tau}; with
Res f4 ~ dphi^{-2} and the celestial kernel ~ dphi^{5/2}, the corner integrand is dphi^{tau - 3/2}: integrable iff tau > 1/2.
Explicit structures: the (13) leg of a twist-1 exchange in light-ray polarisation, D1^m D3^n X^{-1/2} (D = z.d), m, n in {0, 1, 2}
(these are the exact singular factors of the free-field realisation of the twist-1 family, which contains T): closed forms below.
U_mn(dphi) = (pi^2/64) (1 - X^2)/X * S_mn  (edge residue of t^3 t^3 f4 x structure; the smooth legs to 2 and 4 set to 1).
Usage: python3 v4d_og.py -> cache/v4d_og.txt (~3 min)
"""
import os, sys, time
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '2')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import v4d_core as D

CACHE = os.path.join(HERE, 'cache')
L = []


def pr(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); L.append(s)


def S_mn(dph, eps, m, n, sgn=+1):
    c = np.cos(dph); y = 1 - c
    r = -y + 0j
    p = y/2 - 1j*sgn*eps
    q = y/2 - 1j*sgn*eps
    X = -y/2 + 2j*sgn*eps + eps**2
    h = X**-0.5
    tab = {(0, 0): h, (1, 0): -p*h/X, (0, 1): q*h/X, (1, 1): (-3*p*q/X + r)*h/X, (2, 0): 3*p**2*h/X**2, (0, 2): 3*q**2*h/X**2,
           (2, 1): 3*(-2*p*r*X**-2.5 + 5*p**2*q*X**-3.5),
           (1, 2): 3*(2*q*r*X**-2.5 - 5*p*q**2*X**-3.5),
           (2, 2): 3*(2*r**2*X**-2.5 - 20*p*q*r*X**-3.5 + 35*p**2*q**2*X**-4.5)}
    return tab[(m, n)], X


def U_mn(dph, eps, m, n, sgn=+1):
    S, X = S_mn(dph, eps, m, n, sgn)
    return (np.pi**2/64)*(1 - X**2)/X*S


def main():
    t0 = time.time()
    pr('# v4d_og: u-channel twist-1 exchange (light-ray polarised, exact (13)-leg structures) at the corner')
    # ---- (1) scaling checks of the structures
    pr('(1) contact-region scaling: |S_mn(dphi = s sqrt(eps))| ~ eps^P (fit eps = 1e-4..1e-8, s = 1); pointwise |S_mn(dphi, eps=1e-14)| ~ dphi^Q')
    for (m, n) in ((0, 0), (1, 0), (1, 1), (2, 0), (2, 1), (2, 2)):
        e = np.array([1e-4, 1e-6, 1e-8])
        v = np.array([abs(S_mn(np.sqrt(x), x, m, n)[0]) for x in e])
        P = np.polyfit(np.log(e), np.log(v), 1)[0]
        dd = np.array([1e-2, 1e-3, 1e-4])
        w = np.array([abs(S_mn(x, 1e-14, m, n)[0]) for x in dd])
        Q = np.polyfit(np.log(dd), np.log(w), 1)[0]
        pr(f'    (m,n) = ({m},{n}): P = {P:+.3f} (lemma -1/2);  Q = {Q:+.3f} (lemma -1);  identically zero? {np.all(w == 0)}')
    # ---- (2) raw per-nu corner pairing: finite eps -> 0 limit, approach eps^{1/4}
    pr('(2) raw per-nu pairing G_nu(eps) = int C_nu U_mn (exact Fourier pairing): eps -> 0 behaviour')
    EP = (1e-3, 1e-4, 1e-5, 1e-6, 1e-7)
    for (m, n) in ((1, 1), (2, 2), (0, 0)):
        for nu in (0.5, 2.0, 6.0):
            g = []
            for e in EP:
                N = 2**18 if e < 3e-7 else (2**17 if e < 3e-6 else 2**15)
                g.append(D.pair_exact(nu, lambda x: U_mn(x, e, m, n), N=N))
            g = np.array(g)
            dg = np.abs(np.diff(g))
            ex = np.log(dg[1:]/dg[:-1])/np.log(10.0)
            pr(f'    (m,n)=({m},{n}) nu={nu}: |G| = ' + ' '.join(f'{abs(x):.4e}' for x in g) +
               '  ; |G(eps_k) - G(eps_k+1)| ratio exponents = ' + ' '.join(f'{-x:+.3f}' for x in ex) + '  (eps^{1/4}: +0.25)')
    # ---- (3) nu-dependence of the per-nu pairing and heavy normalisation
    pr('(3) per-nu size |G_nu(eps = 1e-5)| and heavy-normalised size (unit: v4c kernel at nu = 0), structure (1,1):')
    lN0 = D.logNcom(0.0)
    for nu in (1, 3, 5, 8, 12, 16, 20, 25):
        v = D.pair_exact(float(nu), lambda x: U_mn(x, 1e-5, 1, 1), N=2**17)
        pr(f'    nu = {nu:3d}: |G| = {abs(v):.3e}   |G| N0/N = {abs(v)*np.exp(lN0 - D.logNcom(nu)):.3e}')
    # ---- (4) de-shifted with the actual weights (XP value / deriv), truncated at nu_c
    pr('(4) de-shifted pairing with the actual weights, sum_{nu <= nu_c} dnu chi N0/N G_nu(eps), structure (1,1):')
    Wt = D.weights()
    for nm in ('xp_value', 'xp_deriv'):
        nu, chi, dnu = Wt[nm]
        w = chi*np.exp(lN0 - D.logNcom(nu))
        sel = nu <= 20.0
        Gs = {}
        for e in (1e-3, 1e-5, 1e-7):
            N = 2**18 if e < 3e-7 else 2**17
            Gs[e] = np.array([D.pair_exact(float(v), lambda x: U_mn(x, e, 1, 1), N=N) for v in nu[sel]])
        for nuc in (4, 8, 12, 16, 20):
            m_ = nu[sel] <= nuc
            vals = [abs(dnu*np.sum(w[sel][m_]*Gs[e][m_])) for e in (1e-3, 1e-5, 1e-7)]
            pr(f'    {nm} nu_c = {nuc:3d}: |G_ds| = ' + '  '.join(f'{v:.3e}' for v in vals) + '   (eps = 1e-3 / 1e-5 / 1e-7)')
    pr(f'[{time.time() - t0:.0f}s]')
    open(os.path.join(CACHE, 'v4d_og.txt'), 'w').write('\n'.join(L) + '\n')


if __name__ == '__main__':
    main()
