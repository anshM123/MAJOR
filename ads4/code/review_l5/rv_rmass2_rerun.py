"""
rt_rmass2.py: per-unit-cell certified bounds of int |r| (r = realised repaired density - psi2) on [940, mu_end], by EPSW_CERT's own panel
machinery (epswcert/ew_cert.py panels/f_series/f_disc_sup/abs_poly_bounds, imported read-only; same N = 40, Q = 5, rho, error terms),
regrouped by unit cells, plus EPSW's Gaussian tail beyond mu_end.  Output cache/rt_rmass2.json: [[T, T+1, mass], ...].
"""
import os, sys, json
EPS = '/home/user/MAJOR/ads4/code/epswcert'
sys.path.insert(0, EPS)
sys.dont_write_bytecode = True
import ew_common as C
import ew_cert as EC
from flint import arb, arb_poly
from fractions import Fraction as Fr
from an_common import CACHE

def main():
    W = C.Weight(1000, os.path.join(C.CACHE, 'k2_continuum_D_1000_rep.json'))
    P, mu_end = EC.panels(W)
    N, Q = EC.N, EC.Q
    cells = {}
    for ip, (a, b, inside) in enumerate(P):
        af = float(a.lower())
        if af < 897.9:
            continue
        m = (a + b)/2; h = arb(((b - a)/2).upper())
        mf = float(m.mid())
        rho = min(0.5, 0.6*mf)
        rad = float(m.rad()) + 1e-30
        Mf, Ff = C.f_disc_sup(W, mf, rho + rad)
        x = h/rho
        fc = C.f_series(W, m, N)
        E = [2*c.real for c in fc]
        if inside:
            ps = W.psi_taylor(m, N)
            rc = [E[n] - ps[n] for n in range(N)]
        else:
            rc = E
        Mn = Mf + Ff
        er = 2*(Mn*x**N/(1 - x) + Ff)
        pr_ = arb_poly(rc)
        hs = h/Q
        i1 = arb(0)
        for j in range(Q):
            cj = -h + (2*j + 1)*hs
            s1, j1 = EC.abs_poly_bounds(EC.coeffs(EC.shift(pr_, cj), N), hs)
            i1 += j1 + 2*hs*er
        T = int(af)                      # panel [a, b] with b - a <= 1/10 lies in [T, T + 1.1): charge it to cell T (cells overlap at most by one panel)
        cells[T] = cells.get(T, 0.0) + float(i1.upper())
    T1, T2 = EC.tail(W, mu_end, 're')
    out = [[T, T + 1.1, v] for T, v in sorted(cells.items())] + [[float(mu_end), 1e9, float(T1.upper())]]
    json.dump(out, open(os.path.join(CACHE, 'rt_rmass2_rerun.json'), 'w'))
    for r in out:
        if r[0] in (940, 980, 990, 995, 998, 999, 1000, 1002, 1005, 1010, 1020, 1028, 1035, 1040, 1058) or r[1] > 1e8:
            print(r)
    print('total', sum(r[2] for r in out))

if __name__ == '__main__':
    main()
