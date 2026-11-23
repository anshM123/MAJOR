"""
p2con/pc_p4.py (P2_CONTRACT.md §2): the flat nu^6 tower of P4 and the hybrid floor.

(a) Flat scaling limit of the AdS Regge-moment series.  In the double limit m -> inf, nu^2/m^2 = y and nu beta = b fixed, a ring
    structure nu^{2j} (eta^2-1)^i P^{(n)} survives at O(y^k) iff 2j - 2i + 2n = 2k.  Using P_nu(cosh beta) -> J0(b),
    (eta^2-1) P' -> -b J1(b), (eta^2-1)^2 P'' reduced by Legendre's equation, CHMRS's Delta_phi-INDEPENDENT top terms give
        O(y):   g1(b) = (11/6) J0 - (1/6) b J1                                        (B, C2 at eta = 1; = the exact EE result, §1)
        O(y^2): g2(b) = (3594/1440) J0 - (11/30) b J1 - (1/72) b^2 J0                  (nu^4 P|_{eta=1}, nu^4 eta(eta^2-1)P', nu^6 (eta^2-1) P)
    These are compared with the y-expansion at fixed b of the FLAT partial-wave image used by the hybrid floor
        h(y; J) = (2 - y) Dtilde_J(y) / 2,   Dtilde = P_J(1-2y) (++),  P^{(0,8)}_{J-4}(1-2y) (+-)   (omegau/nf_common.flat_parts_nf)
    and with the unstripped k = 2 envelope E_2(y) = (2-y)/(1-y)^2 (CHMRS (3.40)): numerically, by extrapolation in J at fixed b.
(b) The beyond-hybrid remainder at O(1/m^4): the hybrid contains the flat image exactly, so P4 has to bound only
        R4 = [Delta_phi-dependent / AdS part of the O(1/m^4) series]  +  [mismatch of the flat tower between AdS (beta, nu, P) and flat
             (bhat, p, J0) kinematics],
    computed here per sample state as exact Lemma-I costs relative to the K4 unit (4a = eta^2 nu^4 P), together with the O(1/m^2)
    analogue (the beyond-hybrid part of the nu^2 classes B, C2: the 'nu^2 x beta^2' scenario of LARGESPIN_W, now exact).
Output: cache/p4.txt, cache/p4.json
"""
import os, sys, json, io
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'lsw')); sys.path.insert(0, os.path.join(HERE, '..', 'omegau'))
import numpy as np
from scipy.special import j0, j1, eval_legendre, eval_jacobi
import lsw_common as W

OUT = io.StringIO()
pr = lambda s: (print(s, flush=True), OUT.write(s + '\n'))


def Dt(h, J, y):
    x = 1 - 2*y
    return eval_legendre(J, x) if h == '++' else eval_jacobi(J - 4, 0, 8, x)


def g1(b):
    return (11/6)*j0(b) - b*j1(b)/6


def g2(b):
    return (3594/1440)*j0(b) - (11/30)*b*j1(b) - b*b*j0(b)/72


def expand_flat(h, b, env):
    """fit h(y)/1 - J0(b) = c1 y + c2 y^2 + c3 y^3 + c4 y^4 at fixed b over J (y = (b/(2J+1))^2)."""
    Js = np.array([60, 80, 110, 150, 200, 280, 400, 560, 800, 1100, 1600, 2300])
    ys = (b/(2*Js + 1.0))**2
    vals = np.array([(2 - y)*env(y)*Dt(h, J, y)/2 - j0(b) for J, y in zip(Js, ys)])
    A = np.vstack([ys**k for k in range(1, 6)]).T
    c, *_ = np.linalg.lstsq(A/ys[:, None], vals/ys, rcond=None)
    return c[:2]


def main():
    pr(__doc__.split('Output')[0])
    pr('(a) y-expansion of the flat image at fixed b = (2J+1) sqrt(y):  [c1, c2] vs AdS top terms [g1, g2]')
    envs = {'stripped (2-y) Dt (hybrid floor)': lambda y: 1.0, 'unstripped E_2 = (2-y)/(1-y)^2': lambda y: (1 - y)**-2}
    res = {}
    for name, env in envs.items():
        for h in ('++', '+-'):
            rows = []
            for b in (0.5, 1.3, 2.7, 4.1, 6.0):
                c1, c2 = expand_flat(h, b, env)
                rows.append((b, c1, g1(b), c2, g2(b)))
            dev1 = max(abs(r[1] - r[2]) for r in rows); dev2 = max(abs(r[3] - r[4]) for r in rows)
            pr(f'   {name:34s} {h}: max|c1 - g1| = {dev1:.2e}, max|c2 - g2| = {dev2:.2e};  e.g. b = 2.7: c1 = {rows[2][1]:+.6f} '
               f'(g1 {rows[2][2]:+.6f}), c2 = {rows[2][3]:+.6f} (g2 {rows[2][4]:+.6f})')
            res[f'{name}|{h}'] = dict(dev1=dev1, dev2=dev2, rows=rows)
    # (b) beyond-hybrid remainders per state
    R = json.load(open(os.path.join(HERE, '..', 'lsw', 'cache', 'costs_1000.json')))
    nu = W.NU
    pr('\n(b) beyond-hybrid parts, exact Lemma-I costs (de-shift form, ring weight 2) per unit m^-4 resp. m^-2, in units of the K4 '
       'representative 4a = eta^2 nu^4 P (cost c4a):')
    pr('   J   m/L   beta | K4-equivalent of: nu^6 tower raw (P2 §7.1) | nu^6 mismatch (AdS-flat) | full O(y^2) flat-tower mismatch '
       '| Delta_phi-dep. AdS nu^4 P (eta^2-1) | sum (remainder, sign-free) | O(1/m^2): beyond-hybrid nu^2 (B+C2) cost / H')
    rows = []
    for r in R:
        b = r['beta']; m = r['m']; mfl = r['mfl']; J = r['J']
        P, dPb = W.P_and_dP(b); P = P.real; dPb = dPb.real
        ch, sh = np.cosh(b), np.sinh(b)
        dPe = dPb/sh
        bh = (2*J + 1)/mfl; p = nu/W.L
        bb = bh*p
        act = lambda f: W.act(f)[0]
        c4a = abs(act(ch**2*nu**4*P))
        nu6 = act(-(ch*ch - 1)*nu**6*P/72)
        mis6 = act(-(nu**4)*((ch*ch - 1)*nu**2*P - bb**2*j0(bb))/72)
        top4 = (3594/1440)*nu**4*P + (11/30)*ch*(ch*ch - 1)*nu**4*dPe - (ch*ch - 1)*nu**6*P/72
        flat4 = nu**4*g2(bb)
        misall = act(top4 - flat4)
        ads4 = act(5189*(ch*ch - 1)*nu**4*P/1440)
        # O(1/m^2): beyond-hybrid nu^2 classes: AdS (11/6 eta nu^2 P + 1/6 (eta^2-1) nu^2 P') - flat nu^2 g1(bhat p)
        mis2 = act((11/6)*ch*nu**2*P + (1/6)*(ch*ch - 1)*nu**2*dPe - nu**2*g1(bb))
        full2 = abs(act((11/6)*ch*nu**2*P)) + abs(act((1/6)*(ch*ch - 1)*nu**2*dPe))
        f = {h: v['floor'] for h, v in r['fl'].items()}
        H = f['++']*f['+-']/(f['++'] + f['+-']) if len(f) == 2 else list(f.values())[0]
        row = dict(J=J, mfl=mfl, beta=b, m=m, H=H, c4a=c4a, K_nu6=abs(nu6)/c4a, K_mis6=abs(mis6)/c4a, K_misall=abs(misall)/c4a,
                   K_ads4=abs(ads4)/c4a, mis2=abs(mis2)/m**2/H, full2=full2/m**2/H)
        row['K_rem'] = row['K_misall'] + row['K_ads4']
        rows.append(row)
        pr(f"  {J:3d} {mfl:.3f} {b:.4f} | {row['K_nu6']:9.3g} | {row['K_mis6']:9.3g} | {row['K_misall']:9.3g} | {row['K_ads4']:9.3g} | "
           f"{row['K_rem']:9.3g} | {row['mis2']:.3g} (full nu^2 classes: {row['full2']:.3g})")
    json.dump(dict(expansion=res, rows=rows), open(os.path.join(HERE, 'cache', 'p4.json'), 'w'))
    open(os.path.join(HERE, 'cache', 'p4.txt'), 'w').write(OUT.getvalue())


if __name__ == '__main__':
    main()
