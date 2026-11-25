"""kc_bounds.py -- certified headline numbers of item (1) as a function of the |g3|^2 light coefficient K (ads4/K_CERT.md).

tau_0(K) = (240/Lam^4) sqrt(E/K),  E = realised Einstein light value of Omega^EE (= 16 Lam V3; equals 16 L_A to 1e-9 for the continuum
weights, 1e-5 for the grids), K = K0 + dK0 + H17.  Everything below is ball arithmetic; printed upper bounds are rounded UP.
Inputs: kc_dk0.main (continuum weights, 300 / 1000) and kc_grid.main (grid weights, 2000 / 5000).
Also: the H17 sensitivity coefficients (moments) and the model band 36/Lam^2 * sum |c| D int psi p^(2D-2) recomputed exactly.
Usage: python3 kc_bounds.py   -> cache/kc_bounds.txt
"""
import sys, os, math
from fractions import Fraction as Fr
from flint import arb, ctx
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kc_common import *  # noqa
import kc_dk0, kc_grid

ctx.prec = 128


def up(x, sig=4):
    """round the upper end of ball x up to `sig` significant digits (decimal string)."""
    u = float(x.upper()) * (1 + 1e-15)
    e = math.floor(math.log10(abs(u))) - sig + 1
    return '%.*e' % (sig - 1, math.ceil(u / 10 ** e) * 10 ** e)


def lo(x, sig=4):
    u = float(x.lower()) * (1 - 1e-15)
    e = math.floor(math.log10(abs(u))) - sig + 1
    return '%.*e' % (sig - 1, math.floor(u / 10 ** e) * 10 ** e)


def main():
    out = []

    def pr(s):
        print(s, flush=True); out.append(s)
    p2, p3 = nf02()
    K0 = arb(fq(K0_exact()))
    # model band (h3_H17_model.py heuristic, recomputed exactly): |dk| = 36/Lam^2 [2*3*int psi2 p^4 + 2*2*int psi3 p^2]
    band_num = 36 * (2 * 3 * poly_moment(p2, 4) + 2 * 2 * poly_moment(p3, 2))
    pr('# kc_bounds: model H17 band on K = %.6f/Lam^2 (bound factor (1 -+ dk)^(-1/2) ~ 1 +- %.4f/Lam^2)' % (float(band_num), float(band_num) / 2))
    pr('# K0 = 1 + %.4e (exact rational)' % float(K0_exact() - 1))
    rows = {}
    for Dg in (300, 1000):
        _, r = kc_dk0.main(Dg)
        rows[Dg] = dict(kind='continuum', dK0=r['dK0'], dK0s=r['dK0s'], E=16 * r['Lam'] * r['V3'], LA=r['LA'], Lam=r['Lam'],
                        M=[r['M0'], r['M1'], r['M2'], r['M3']])
    _, g = kc_grid.main([2000, 5000])
    for Dg in (2000, 5000):
        r = g[Dg]
        rows[Dg] = dict(kind='grid', dK0=r['dK0'], dK0s=r['dK0s'], E=16 * r['Lam'] * r['V3'], LA=r['LA'], Lam=r['Lam'], M=None)
    targets = {300: ['3.135e-6', '3.14e-6', '3.15e-6', '3.2e-6'], 1000: ['2.892e-8', '2.90e-8', '2.95e-8', '3.0e-8'],
               2000: ['1.929e-9', '1.93e-9', '1.95e-9'], 5000: ['5.327e-11', '5.33e-11', '5.4e-11']}
    pr('')
    pr('| Dg | weight | Lam^2 dK0 (h1) | Lam^2 dK0 (C\') | K_c = K0+dK0 | E = realised 16 L_A | tau0(K=1) | tau0(K_c) | tau0 over model band | ')
    for Dg, r in rows.items():
        L = r['Lam']; L2 = L * L
        tau = lambda K: 240 / L ** 4 * (r['E'] / K).sqrt()
        Kc = K0 + r['dK0']
        Kmin = min_arb(K0 + r['dK0'], K0 + r['dK0s']) - arb(fq(band_num)) / L2
        Kmax = max_arb(K0 + r['dK0'], K0 + r['dK0s']) + arb(fq(band_num)) / L2
        pr('| %d | %s | %s | %s | %s | %s | %s | %s | [%s, %s] |' % (
            Dg, r['kind'], (r['dK0'] * L2).str(7, radius=False), (r['dK0s'] * L2).str(3, radius=False), Kc.str(12, radius=False),
            r['E'].str(12, radius=False), tau(K0 * 0 + 1).str(9, radius=False), tau(Kc).str(9, radius=False), lo(tau(Kmax), 5), up(tau(Kmin), 5)))
        r.update(Kc=Kc, Kmin=Kmin, Kmax=Kmax, tau=tau)
    pr('')
    pr('Tolerance of rounded headlines to H17: tau0 <= T  <=>  K >= K*(T) = (240/Lam^4)^2 E/T^2;  h* := K* - K_c (H17 needed), in units of')
    pr('the model band b = %.2f/Lam^2 and of 1/Lam^2.' % float(band_num))
    for Dg, r in rows.items():
        L = r['Lam']; L2 = L * L
        for T in targets[Dg]:
            Tq = arb(T)
            Ks = (240 / L ** 4) ** 2 * r['E'] / (Tq * Tq)
            h = Ks - r['Kc']
            pr('  Dg = %d: tau0 <= %s  <=>  K >= %s ; holds iff H17 >= %s = %s/Lam^2 = %s x model band' % (
                Dg, T, Ks.str(10, radius=False), h.str(6, radius=False), (h * L2).str(6, radius=False), (h * L2 / arb(fq(band_num))).str(4, radius=False)))
    pr('')
    pr('H17 sensitivity (continuum weights). If -C\'_EE|_{g3} = 2piG|g3|^2 (nu/Lam)^6 [1 + a4/nu^2 + a2/nu^4 + a0/nu^6] (k = 2, h1 reading)')
    pr('and -C_3|_{g3} = -2piG|g3|^2 (nu/Lam)^4 [1 + b2/nu^2 + b0/nu^4], then exactly')
    pr('   H17 = -(2/Lam^2)[a4 M2 + a2 M1/Lam^2 + a0 M0/Lam^4] + (2/Lam^2)[b2 N1 + b0 N0/Lam^2],  M_j = int Re X_theta p^(2j) dnu/Lam,  N_j = int psi3 p^(2j)')
    for j in range(3):
        pr('   N_%d = %.15f (exact rational)' % (j, float(poly_moment(p3, 2 * j))))
    for Dg in (300, 1000):
        M = rows[Dg]['M']
        pr('   Dg = %d: M0 = %s, M1 = %s, M2 = %s, M3 = %s' % (Dg, M[0].str(12), M[1].str(12), M[2].str(12), M[3].str(12)))
    pr('   leading:  Lam^2 H17 = -2 M2 a4 + 2 N1 b2 + O(Lam^-2) = -%.4f a4 + %.4f b2' % (2 * float(rows[1000]['M'][2]), 2 * float(poly_moment(p3, 2))))
    pr('')
    pr('Quartic upper end: K does not enter (2a); in (2b) it enters only through the sign condition kappa + eps K/8 > 0 (c1 = 0).')
    return out


def min_arb(a, b):
    return a.min(b)


def max_arb(a, b):
    return a.max(b)


if __name__ == '__main__':
    o = main()
    open(os.path.join(HERE, 'cache', 'kc_bounds.txt'), 'w').write('\n'.join(o) + '\n')
