"""agproof/ag_edge2.py (AG_PROOF.md §5.3): the exact scalar's (EDGE-R) quantity E_chi[G] = int_{p_w}^1 psi_2 chi G dp (k = 2), with the
SAME cutoffs (p_e, sigma = s T, T = t 25/(0.97 Lam beta)) as Lemma P4W-chi's grid, at the exact gap states J = 40, 48, 64, 100, in units
of S_E2 := 2 V_Fa2(1) e_sc(1) (the k = 2 edge cost of a flat-profile wave pair of amplitude e_sc(1), Lemma F).  The quantity includes
the part of the cutoff region that overlaps the window (chi(p_w) > 0); it is exactly what Lemma P4W-chi's premise (EDGE-R) bounds.
-> cache/edge2.txt, cache/edge2.json
"""
import os, sys, json
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.special import erfc, j0, eval_legendre
import ag_exact as AE
import ag_flat as AF
import p4w_common as C          # read-only
import p4w_contour as K         # read-only

L = AE.LAM
OUT = []
pr = lambda s: (print(s, flush=True), OUT.append(s))
GRID_PE = (0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75)
GRID_T = (1.0, 0.8, 0.65, 0.5)
GRID_S = (0.35, 0.45, 0.6)


def main():
    pr(__doc__.split('->')[0])
    res = {}
    for J in (40, 48, 64, 100):
        st = AE.state(J); m = st['m']; beta = st['beta']; mfl = st['mfl']
        L3 = AE.build(st, numax=L*1.02, prec=400, conf=(320, 960, 0.5))
        k0 = 0.97*beta*L; T0 = 25.0/k0
        pw = 3.0/beta/L
        kmax = (2*J + 1)/mfl/np.sqrt(1 - (1/mfl)**2)
        nseg = int(np.ceil((1 - pw)*1.2*kmax/np.pi*2)) + 40
        xs, ws = np.polynomial.legendre.leggauss(16)
        br = np.linspace(pw, 1.0, nseg + 1)
        p = np.concatenate([a + (b - a)*(xs + 1)/2 for a, b in zip(br[:-1], br[1:])])
        wp = np.concatenate([ws*(b - a)/2 for a, b in zip(br[:-1], br[1:])])
        nu = L*p; x = nu/m
        ex = np.array([float(L3.value_fast(2, float(v)).mid()) for v in nu])/2
        flat = (2 - x*x)/(1 - x*x)**2*eval_legendre(J, 1 - 2*x*x)/2
        cur = C.Pd(beta, nu)[0].real - j0((2*J + 1)/mfl*p)
        G = ex - flat - cur
        ps = K.psi_c(2, p).real
        fl = AF.edge_compare(J, mfl)
        esc1 = float(K.e_sc(beta, mfl)(np.array([1.0]))[0])
        SE2 = 2*fl['VF2']*esc1
        kl_e = lambda pe: (2*J + 1)/mfl/np.sqrt(1 - (pe/mfl)**2)
        rows = []
        for pe in GRID_PE:
            for tf in GRID_T:
                for sf in GRID_S:
                    sg = sf*tf*T0
                    if kl_e(pe)*sg < 10:
                        continue
                    chi = erfc((pe - p)/sg)/2
                    E = float(np.sum(wp*ps*chi*G))
                    rows.append(dict(pe=pe, Tfrac=tf, s=sf, sigma=sg, E=E, ratio=abs(E)/SE2, chi_pw=float(erfc((pe - pw)/sg)/2)))
        res[J] = dict(J=J, mfl=mfl, SE2=SE2, rows=rows, total_complement=float(np.sum(wp*ps*G)))
        rr = np.array([r['ratio'] for r in rows])
        pr(f"J = {J}: S_E2 = {SE2:.3e}; |E_chi[G]|/S_E2 over the grid ({len(rows)} cutoffs): min {rr.min():.3f}, median {np.median(rr):.3f}, "
           f"max {rr.max():.3f}; full complement int_(p_w)^1 psi2 G = {res[J]['total_complement']:+.3e} (= {abs(res[J]['total_complement'])/SE2:.3f} S_E2)")
        for r in sorted(rows, key=lambda r: -r['ratio'])[:3]:
            pr(f"     largest: p_e = {r['pe']}, T/T0 = {r['Tfrac']}, sigma = {r['sigma']:.3f} (chi(p_w) = {r['chi_pw']:.0e}): {r['ratio']:.3f}")
    json.dump({str(k): v for k, v in res.items()}, open(os.path.join(AE.CACHE, 'edge2.json'), 'w'), default=float)
    open(os.path.join(AE.CACHE, 'edge2.txt'), 'w').write('\n'.join(OUT) + '\n')


if __name__ == '__main__':
    main()
