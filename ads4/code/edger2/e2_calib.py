"""
edger2/e2_calib.py (EDGE_R2.md §4): calibration of the remainder treatment on the EXACT scalar (Delta_phi = 3, k = 2; agproof/ag_exact
Arb values, read-only; the er_scalar.py conventions: G = ex/2 - E_2 P_J/2 - (P - J0), Gamma = m^2 G/(E_2/2), weight psi_2 chi dp).
For each even J of the batch-3 states (twist edge of Delta = 1000 + J = the W state), with the GRAVITON's own node layout
(er_jobs.state_nodes: bulk 0.025 grid, two GK15 windows [a, x_mid], [x_mid, x1]) and the graviton's last measured n1 node x_l
(e2_merge: nu/Delta < 0.93):
  (1) E_full  = E_chi[G] (dense quadrature),
      E_det   = the analogue of what the n1 rows determine: dense integral on [p_w, a) + K15 sum over the measured GK nodes (x <= x_l),
      E_rem   = E_full - E_det (the part carried by the missing GK nodes + the K15 error of the upper window),
      quad_up = K15 sum with ALL upper-window nodes - dense integral over the upper window (GK15 error there);
  (2) the growth envelope used for the graviton's missing nodes, G_env = A ((1 - x_l)/(1 - x))^2 (A = max|Gamma| of the last three
      measured nodes): c_pt = max over the missing GK nodes of |Gamma(x)|/G_env(x)  (the pointwise factor the envelope needs),
      c_int = sum |w||Gamma| / sum |w| G_env  and  c_rem = |E_rem| / sum |w| G_env  (signed remainder vs the sign-free envelope);
  (3) the exponent of the measured |Gamma| growth on [0.8, x_l] (fit of log max-envelope vs log(1 - x)) for comparison with the
      graviton's (e2_report).
Units: S_E2 = 2 V_Fa2(1) e_sc(1) (er_scalar).  -> cache/calib.json / .txt
Usage: python3 e2_calib.py 26,30,34,38,42,46,50
"""
import e2_guard  # noqa: F401
import os, sys, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, os.path.join(CODE, 'edger')); sys.path.insert(0, os.path.join(CODE, 'highspin'))
import numpy as np
from scipy.special import j0, eval_legendre
import er_core as R
import ag_exact as AE
import ag_flat as AF
import p4w_common as C
import p4w_contour as K
import hs_core as HC
import er_jobs as EJ
CACHE = os.path.join(HERE, 'cache')
L = AE.LAM
XL = json.load(open(os.path.join(CACHE, 'xl.json'))) if os.path.exists(os.path.join(CACHE, 'xl.json')) else {}


def main(Js):
    t0 = time.time(); OUT = []
    pr = lambda s='': (print(s, flush=True), OUT.append(s))
    pr(__doc__.split('Units')[0])
    res = {}
    for J in Js:
        st = AE.state(J); m = st['m']; beta = st['beta']; mfl = st['mfl']
        L3 = AE.build(st, numax=L*1.02, prec=400, conf=(320, 960, 0.5))
        ex = lambda nu: np.array([float(L3.value_fast(2, float(v)).mid()) for v in np.atleast_1d(nu)])/2

        def Gof(x):
            x = np.atleast_1d(np.asarray(x, float)); nu = x*m; p = nu/L
            E2h = (2 - x*x)/(1 - x*x)**2/2
            return ex(nu) - E2h*eval_legendre(J, 1 - 2*x*x) - (C.Pd(beta, nu)[0].real - j0((2*J + 1)/mfl*p)), E2h
        wr = R.w_state(J, mfl)
        s = EJ.state_nodes(J, wr['mfl'], wr['pe'], wr['sigma'], wr['pw'])
        x1 = s['x1']; a = 0.8
        for _ in range(4):                       # the windows actually run (e2_gamma.layout: er_jobs reports edge_a one shift late)
            xm = 0.5*(a + x1); e = EJ.gk(a, xm) + EJ.gk(xm, x1); a_used = a
            if min(abs(np.array(e) - s['x0'])) >= 0.004:
                break
            a += 0.006
        a = a_used; xm = 0.5*(a + x1)
        pe, sg = wr['pe'], wr['sigma']
        pw = 3.0/beta/L
        xl = XL.get(str(J), max(x for x in s['n1']))
        # dense quadrature on [p_w, 1] (16-pt GL panels)
        kmax = (2*J + 1)/mfl/np.sqrt(1 - (1/mfl)**2)
        nseg = int(np.ceil((1 - pw)*1.2*kmax/np.pi*2)) + 40
        gx, gw = np.polynomial.legendre.leggauss(16)
        br = np.unique(np.concatenate([np.linspace(pw, 1.0, nseg + 1), [a*mfl, xm*mfl]]))   # panel edges at the window ends
        p = np.concatenate([u + (v - u)*(gx + 1)/2 for u, v in zip(br[:-1], br[1:])])
        wp = np.concatenate([gw*(v - u)/2 for u, v in zip(br[:-1], br[1:])])
        x = p*L/m
        G, E2h = Gof(x)
        ps = K.psi_c(2, p).real
        ch = R.chi(p, pe, sg)
        fl = AF.edge_compare(J, mfl)
        esc1 = float(K.e_sc(beta, mfl)(np.array([1.0]))[0])
        SE2 = 2*fl['VF2']*esc1
        f = wp*ps*ch*G
        E_full = float(f.sum())
        E_below_a = float(f[x < a - 1e-12].sum())
        up = (x >= xm)
        E_up_dense = float(f[up].sum())
        # GK windows (weights in x -> dp = mfl dx)
        dets, remw, allup, alllo = 0.0, [], 0.0, 0.0
        Gn = {}
        for (lo, hi) in ((a, xm), (xm, x1)):
            xe, wk, wg = HC.gk15_nodes(lo, hi)
            Ge, E2e = Gof(xe)
            pe_ = xe*mfl
            we = wk*mfl*K.psi_c(2, np.minimum(pe_, 1)).real*R.chi(pe_, pe, sg)*(pe_ <= 1)
            for xi, wi, gi, e2i in zip(xe, we, Ge, E2e):
                Gn[float(xi)] = (float(gi), float(e2i))
                if xi <= xl + 1e-9:
                    dets += wi*gi
                else:
                    remw.append((float(xi), float(wi), float(gi), float(e2i)))
                if lo == xm:
                    allup += wi*gi
                else:
                    alllo += wi*gi
        E_det = E_below_a + dets
        E_rem = E_full - E_det
        quad_up = allup - E_up_dense
        quad_lo = alllo - float(f[(x >= a) & (x < xm)].sum())
        # envelope: Gamma = m^2 G / E2h at the last three measured GK / grid nodes (<= x_l)
        cand = sorted([xx for xx in Gn if xx <= xl + 1e-9] + [float(v) for v in s['n1'] if v >= 0.8 and v <= xl + 1e-9])
        last = sorted(set(round(v, 7) for v in cand))[-3:]
        GL, E2L = Gof(np.array(last))
        A = float(np.max(np.abs(m*m*GL/E2L)))
        cpt, num, den = 0.0, 0.0, 0.0
        for (xi, wi, gi, e2i) in remw:
            Gam = m*m*gi/e2i
            env = A*max(1.0, ((1 - xl)/(1 - xi))**2)
            cpt = max(cpt, abs(Gam)/env)
            num += abs(wi)*abs(gi); den += abs(wi)*env*e2i/(m*m)
        # growth exponent of the measured envelope on [0.8, x_l]
        xg = np.linspace(0.8, xl, 400); Gg, E2g = Gof(xg); Gam_g = np.abs(m*m*Gg/E2g)
        nb = 8; edges = np.linspace(0.8, xl, nb + 1)
        mx = [Gam_g[(xg >= u) & (xg <= v)].max() for u, v in zip(edges[:-1], edges[1:])]
        mid = 0.5*(edges[:-1] + edges[1:])
        expo = float(-np.polyfit(np.log(1 - mid), np.log(mx), 1)[0])
        # beyond x_l: true growth exponent of the envelope of |Gamma| on [x_l, x1]
        xg2 = np.linspace(xl, x1 - 1e-4, 400); Gg2, E2g2 = Gof(xg2); Gam2 = np.abs(m*m*Gg2/E2g2)
        e2 = np.linspace(xl, x1 - 1e-4, 6); mx2 = [Gam2[(xg2 >= u) & (xg2 <= v)].max() for u, v in zip(e2[:-1], e2[1:])]
        expo2 = float(-np.polyfit(np.log(1 - 0.5*(e2[:-1] + e2[1:])), np.log(mx2), 1)[0])
        r = dict(J=J, mfl=mfl, x1=x1, xl=xl, a=a, xm=xm, SE2=SE2, E_full=E_full/SE2, E_det=E_det/SE2, E_rem=E_rem/SE2,
                 quad_up=quad_up/SE2, quad_lo=quad_lo/SE2, env_sf=den/SE2, true_sf=num/SE2, c_pt=cpt, c_int=num/den if den else 0.0,
                 c_rem=abs(E_rem)/den if den else 0.0, A=A, nmiss=len(remw), expo_meas=expo, expo_beyond=expo2)
        res[J] = r
        pr(f"J = {J} (m/L = {mfl:.4f}, x1 = {x1:.4f}, x_l = {xl:.4f}, {len(remw)} missing GK nodes): E_full = {r['E_full']:+.3f}, "
           f"E_det = {r['E_det']:+.3f}, E_rem = {r['E_rem']:+.3f} S_E2 (GK15 error: lower window {r['quad_lo']:+.4f}, upper {r['quad_up']:+.4f}); "
           f"envelope A((1-x_l)/(1-x))^2 sign-free {r['env_sf']:.3f}, true sign-free {r['true_sf']:.3f} S_E2; c_pt = {cpt:.2f}, "
           f"c_int = {r['c_int']:.2f}, c_rem = |E_rem|/env = {r['c_rem']:.2f}; |Gamma| growth exponent on [0.8, x_l] {expo:.2f}, "
           f"beyond x_l {expo2:.2f}  [{time.time() - t0:.0f}s]")
        json.dump({str(k): v for k, v in res.items()}, open(os.path.join(CACHE, 'calib.json'), 'w'), default=float, indent=1)
    open(os.path.join(CACHE, 'calib.txt'), 'w').write('\n'.join(OUT) + '\n')


if __name__ == '__main__':
    main([int(v) for v in sys.argv[1].split(',')])
