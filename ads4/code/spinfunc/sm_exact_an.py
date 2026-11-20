"""
sm_exact_an (SMOOTH_BOUND.md §2): analysis of the sm_exact records (exact scalar remainder at near-gap states).
Per state and helicity: rem/N_s per k and per x-band (where the remainder sits: small nu / bulk / light-cone edge), the uniform-c2
equivalent C_eq = rem/sens (sens = sum p2 env2 fac/m^2 per c2 unit), the 'smooth projection' of m^2 rho_2 (least squares of
m^2 rho_2 on smooth functions 1, x^2, x^4 with the chain weight, after projecting out the oscillating span {x^k D, x^k x D'/(2J+1)},
k <= 3, on x <= 0.85) and its integrated effect, and the edge profile m^2 rho_k (1 - x^2)^q fits near x_max.
Usage: python3 sm_exact_an.py [Dg ...]  -> stdout (cache/sm_exact_an.txt via tee)
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
EX = os.path.join(HERE, '..', 'exactskel'); sys.path.insert(0, EX)
import numpy as np
from h13_common import F, SIGMA, Lam_of
from nf_common import envk
from scipy.special import eval_legendre
import sm_exact as S

BANDS = ((0.0, 0.1), (0.1, 0.4), (0.4, 0.7), (0.7, 0.85), (0.85, 1.01))


def Dsh(h, J, x):
    from sm_curvfix import Dsh as D
    return D(h, J, x)


def analyse(rec):
    Dg, J, st = rec['Dg'], rec['J'], rec['st']
    mfl, m = st['mfl'], st['m']; u = 1/mfl**2
    n = rec['npn']; pn, pw = rec['pn'], rec['pw']; x = pn/mfl; y = x*x
    r2, r3 = rec['rho'][2][:n], rec['rho'][3][:n]
    p2 = F.psi2(pn)*pw; p3 = F.psi3(pn)*pw
    out = []
    for h, o in rec['hel'].items():
        s = SIGMA[h]; fac = (1 - y)**2 if h == '++' else (1 - y)**-2
        w2 = p2*envk(2, x)*fac; w3 = p3*u*s*envk(3, x)*fac
        Ns = o['Ns']
        band = []
        for a, b in BANDS:
            sel = (x >= a) & (x < b)
            band.append(((w2*r2)[sel].sum()/Ns, (w3*r3)[sel].sum()/Ns))
        # smooth projection of G = m^2 rho_2 on x <= 0.85
        G = m*m*r2; sel = x <= 0.85
        D = Dsh(h, J, x); hh = 1e-6
        Dp = (Dsh(h, J, x + hh) - Dsh(h, J, x - hh))/(2*hh)*x/(2*J + 1)
        osc = [x**k*D for k in range(4)] + [x**k*Dp for k in range(4)]
        smo = [np.ones_like(x), y, y*y]
        A = np.array(osc + smo).T
        wt = np.sqrt(np.abs(w2)) + 1e-300
        c, *_ = np.linalg.lstsq((A*wt[:, None])[sel], (G*wt)[sel], rcond=None)
        rs = np.array(smo).T @ c[8:]
        smooth_eff = (w2*rs/m**2)[sel].sum()/Ns
        out.append((h, Ns, o, band, c[8:], smooth_eff))
    return out


if __name__ == '__main__':
    Dgs = [int(a) for a in sys.argv[1:]] or [300, 1000]
    for Dg in Dgs:
        pk = os.path.join(HERE, 'cache', 'sm_exact_%d.pkl' % Dg)
        if not os.path.exists(pk):
            continue
        db = pickle.load(open(pk, 'rb'))
        for key in sorted(db, key=lambda k: (k[0], str(k[1:]))):
            rec = db[key]; st = rec['st']
            print(f"Dg={Dg} J={rec['J']:3d} key {key} m/L={st['mfl']:.4f} bhat={(2*rec['J']+1)/st['mfl']:.1f} beta={st['beta']:.4f} "
                  f"x_max={1/st['mfl']:.3f}")
            for h, Ns, o, band, cs, se in analyse(rec):
                tot = (o['rem2'] + o['rem3'])/Ns
                print(f"   {h}: rem/N_s {tot:+.2e} [k2 {o['rem2']/Ns:+.2e}, k3 {o['rem3']/Ns:+.2e}]  C_eq {o['Ceq']:+.3f} c2  "
                      f"n-sum {o['nsum_err']/Ns:.0e} odd-k {o['oddk_err']/Ns:.0e}")
                print("        bands (k2/k3, N_s): " + "  ".join(f"[{a:.2f},{b:.2f}) {b2:+.1e}/{b3:+.1e}" for (a, b), (b2, b3) in zip(BANDS, band)))
                print(f"        smooth projection of m^2 rho_2 (x<=0.85): {cs[0]:+.2f} {cs[1]:+.2f} x^2 {cs[2]:+.2f} x^4 c2;  its effect {se:+.2e} N_s")
