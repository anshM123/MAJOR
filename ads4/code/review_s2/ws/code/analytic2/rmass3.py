"""rmass3.py lo hi: per-unit-cell certified bounds of int|r| on [lo, hi) (lo >= 40), by EPSW_CERT's own panel machinery -- an exact copy
of the per-panel computation in code/analytic/rt_rmass2.py (itself EPSW's ew_cert method; imported read-only), only the frequency range
differs.  Output cache/rmass3_lo_hi.json: [[T, T+1.1, mass], ...] (a panel [x, x+1/10] with x in [T, T+1) is charged to cell T)."""
import os, sys, json, time
EPS = '/home/user/MAJOR/ads4/code/epswcert'
sys.path.insert(0, EPS)
sys.dont_write_bytecode = True
import ew_common as C
import ew_cert as EC
from flint import arb, arb_poly
HERE = os.path.dirname(os.path.abspath(__file__))

def main(lo, hi):
    W = C.Weight(1000, os.path.join(C.CACHE, 'k2_continuum_D_1000_rep.json'))
    P, mu_end = EC.panels(W)
    N, Q = EC.N, EC.Q
    cells = {}
    t0 = time.time()
    for ip, (a, b, inside) in enumerate(P):
        af = float(a.lower())
        if af < lo - 1e-9 or af >= hi - 1e-9:
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
        T = int(af)
        cells[T] = cells.get(T, 0.0) + float(i1.upper())
    out = [[T, T + 1.1, v] for T, v in sorted(cells.items())]
    json.dump(out, open(os.path.join(HERE, 'cache', f'rmass3_{lo}_{hi}.json'), 'w'))
    print(f'{len(out)} cells, total {sum(r[2] for r in out):.4e}, [{time.time()-t0:.0f}s]')
    for r in out[::max(1, len(out)//20)]:
        print(r)

if __name__ == '__main__':
    main(int(sys.argv[1]), int(sys.argv[2]))
