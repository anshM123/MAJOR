"""agproof/ag_chi.py (AG_PROOF.md §5): Lemma P4W-chi — the edge handled on the REAL AXIS, the complement with the UNIFORM rate only.

Split the complement cost with the entire cutoff chi(p) = erfc((p_e - p)/sigma)/2:
    2 int_{p_w}^1 psi_2 G = 2 int_0^1 psi_2 chi G  [=: 2 E_chi[G], premise (EDGE-R), real axis]
                           - 2 int_0^{p_w} psi_2 chi G  [<= 2 chi(p_w) int_0^{p_w}|psi_2| sup|G|: included in the window sup bound]
                           + 2 int_{p_w}^1 psi_2 (1 - chi) G  [Cauchy on Q_pm with the entire weight psi_2 (1 - chi)].
For the last term only (AG_U) is used: G = G_+ + G_-, analytic on Q_pm, |G_pm(p +- iy)| <= M(p) e^{-0.97 Lam beta y} (no local rate).
|erfc(a + ib)| <= e^{b^2} erfc(a) <= e^{b^2 - a^2} (a >= 0) bounds the weight on the three sides:
    right side p = 1:  |1 - chi| <= e^{y^2/sigma^2} erfc((1 - p_e)/sigma)/2,
    left side p_w and top: |1 - chi| <= 1 + e^{y^2/sigma^2} erfc(|p_e - p|/sigma)/2 (p < p_e), e^{T^2/sigma^2} erfc((p - p_e)/sigma)/2 (p > p_e).
The edge amplitude enters only through the right-side term, suppressed by erfc((1-p_e)/sigma); the local rate is not used anywhere.
Per state (Theorem W's 47, Theorem L''s), (p_e, sigma, T) are chosen on a grid to minimise the Lam_AG = 1 total; reported:
  * d_chi/H: window + boundary + (1-chi) right side + top at Lam_AG = 1 (plain / +- transferred envelope), i.e. everything but E_chi;
  * eps_E tol = theta H - d_chi, the tolerated |cost of E_chi| (k = 2 + 3), in units of S_E = the same cost for a FLAT-PROFILE wave of
    amplitude e_sc(1) (Lemma F, amplitude included), plain / transferred;
  * the exact scalar's measured |E_chi[G]| in the same units (ag_edge, J = 40, 48, 64, 100).
-> cache/chi.txt, cache/chi.json
"""
import os, sys, json, time
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.special import erfc
import ag_flat as AF
import p4w_common as C          # read-only
import p4w_contour as K         # read-only
import p4w_states as PS         # read-only

L = C.L
_orig = K.contour_bound
PAR = {}
OUT = []
pr = lambda s: (print(s, flush=True), OUT.append(s))


def _V(k, p0, kap, T, wfun, n=200):
    x, w = np.polynomial.legendre.leggauss(n)
    y = (x + 1)/2*T; w = w*T/2
    return float((np.abs(K.psi_c(k, p0 + 1j*y))*np.exp(-kap*y)*wfun(y)*w).sum())


def contour_chi(beta, mfl, pw, Mw, Mfun, mode='uniform', k=2, u=1.0, Aunif=0.0):
    pe, sg, T = PAR['pe'], PAR['sigma'], PAR['T']
    w = 2.0 if k == 2 else u
    k0 = 0.97*beta*L
    left = lambda y: 1 + np.exp(y*y/sg**2)*erfc((pe - pw)/sg)/2
    right = lambda y: np.exp(y*y/sg**2)*erfc((1 - pe)/sg)/2
    bd = 2*w*_V(k, pw, k0, T, left)*Mw
    M1 = float(Mfun(np.array([1.0]))[0])
    ed = 2*w*_V(k, 1.0, k0, T, right)*(M1 + Aunif) + 2*w*_V(k, 1.0, k0, T, right)*Aunif
    x, wq = np.polynomial.legendre.leggauss(500); p = pw + (1 - pw)*(x + 1)/2; wq = wq*(1 - pw)/2
    wt = np.where(p < pe, 1 + np.exp(T*T/sg**2)*erfc((pe - p)/sg)/2, np.exp(T*T/sg**2)*erfc((p - pe)/sg)/2)
    tp = 2*w*float((np.abs(K.psi_c(k, p + 1j*T))*Mfun(p)*wt*wq).sum())*np.exp(-k0*T)
    return bd, ed, tp


_EC = {}


def scale_SE(J, mfl, u, hel):
    """cost (k = 2 + 3) of a flat-profile wave pair of amplitude e_sc(1) (x (1-x1^2)^hel) at the edge: 2*2 V_Fa2 + u*2 V_Fa3."""
    if (J, mfl) not in _EC:
        _EC[(J, mfl)] = AF.edge_compare(J, mfl, ny=0)
    o = _EC[(J, mfl)]
    beta = C.beta_of(J, mfl)
    e1 = float(K.e_sc(beta, mfl, hel=hel)(np.array([1.0]))[0])
    return (2*2*o['VF2'] + u*2*o['VF3'])*e1, o


GRID_PE = (0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75)
GRID_T = (1.0, 0.8, 0.65, 0.5)
GRID_S = (0.35, 0.45, 0.6)          # sigma = s * T


def best_for_state(s, runs, ksh, theta=1.0):
    """one analyse() call records the contour inputs (mode 'local'); the chi-bounds are then evaluated for every grid point and
    combined exactly as p4w_states.analyse combines them (fixed = window + max_e bd; d = fixed + max_e (edge + top))."""
    calls = []

    def rec(beta, mfl, pw, Mw, Mfun, mode='uniform', k=2, u=1.0, Aunif=0.0):
        if mode == 'local':
            calls.append((beta, mfl, pw, Mw, Mfun, k, u, Aunif))
        return 0.0, 0.0, 0.0
    K.contour_bound = rec
    o = PS.analyse(s, runs, ksh)
    K.contour_bound = _orig
    win = o['local']['win']; H = o['H']
    ents = list(o['ent'].keys())
    assert len(calls) == 6*len(ents)
    J, mfl, beta = s['J'], s['mfl'], s['beta']
    T0 = 25.0/(0.97*beta*L)
    kl_e = lambda pe: (2*J + 1)/mfl/np.sqrt(1 - (pe/mfl)**2)
    best = None
    for pe in GRID_PE:
        for tf in GRID_T:
            for sf in GRID_S:
                T = tf*T0; sg = sf*T
                if kl_e(pe)*sg < 10:
                    continue
                PAR.update(pe=pe, sigma=sg, T=T)
                bds, ag_tr, ag_pl = [], [], []
                for i, ek in enumerate(ents):
                    c = calls[6*i:6*i + 6]
                    r = [contour_chi(*a[:5], 'local', k=a[5], u=a[6], Aunif=a[7]) for a in c]
                    # c[0], c[1]: zero envelope k = 2, 3;  c[2], c[3]: transferred k = 2, 3;  c[4], c[5]: plain k = 2, 3
                    bds += [r[2][0] + r[3][0], r[4][0] + r[5][0]]
                    ag_tr.append(r[2][1] + r[2][2] + r[3][1] + r[3][2])
                    ag_pl.append(r[4][1] + r[4][2] + r[5][1] + r[5][2])
                fixed = win + max(bds)
                cand = dict(pe=pe, sigma=sg, T=T, Tfrac=tf, d_plain=(fixed + max(ag_pl))/H, d_tr=(fixed + max(ag_tr))/H, fixed=fixed/H,
                            H=H, u=1/mfl**2, kls=kl_e(pe)*sg, chi_pw=float(erfc((pe - o['pw'])/sg)/2), pw=o['pw'])
                if best is None or cand['d_tr'] < best['d_tr']:
                    best = cand
    return best


def main():
    t0 = time.time()
    pr(__doc__.split('->')[0])
    runs = PS.PJ.load_runs(); ksh = PS.PJ.load_runs('ksh_*_*.json')
    states = json.load(open(os.path.join(C.CODE, 'p2con', 'cache', 'jstar_bh.json')))['states']
    states = sorted(states, key=lambda z: (z['J'] % 2, z['J'], z['mfl']))
    only = os.environ.get('AG_CHI_ONLY')
    if only:
        states = [s for s in states if s['J'] in [int(v) for v in only.split(',')]]
    rows = []
    pr(' J    m/L   | p_e  sigma  T/T0 (kappa_e sigma) | chi(p_w) | d_chi/H plain / tr (Lam_AG = 1) | eps_E tol / S_E: plain / tr')
    for s in states:
        b = best_for_state(s, runs, ksh)
        SEp, o = scale_SE(s['J'], s['mfl'], b['u'], 0.0)
        SEt, _ = scale_SE(s['J'], s['mfl'], b['u'], -2.0)
        tolp = (b['H'] - b['d_plain']*b['H'])/SEp; tolt = (b['H'] - b['d_tr']*b['H'])/SEt
        b.update(J=s['J'], mfl=s['mfl'], SEp=SEp, SEt=SEt, tolE_plain=tolp, tolE_tr=tolt)
        rows.append(b)
        pr(f" {s['J']:3d} {s['mfl']:.3f} | {b['pe']:.2f} {b['sigma']:.3f} {b['Tfrac']:.2f} ({b['kls']:.0f}) | {b['chi_pw']:.0e} | "
           f"{b['d_plain']:.3f} / {b['d_tr']:.3f} | {tolp:9.3g} / {tolt:9.3g}")
    pr('summary: max d_chi/H plain %.3f, transferred %.3f; min eps_E tol/S_E plain %.3g, transferred %.3g' % (
        max(r['d_plain'] for r in rows), max(r['d_tr'] for r in rows), min(r['tolE_plain'] for r in rows), min(r['tolE_tr'] for r in rows)))
    # exact scalar: measured |E_chi[G]| / (2 V_Fa2 e_sc(1)) at the 4 exact states (k = 2 only; ag_edge)
    try:
        ed = json.load(open(os.path.join(AF.CACHE, 'edge.json')))
        for J, o in ed.items():
            Es = [abs(r['EG']) for r in o['rows']]
            fl = AF.edge_compare(int(J), o['mfl'])
            pr(f"  exact scalar J = {J}: |E_chi[G]| = {min(Es):.2e}..{max(Es):.2e}; in units of 2 V_Fa2(1) e_sc(1): "
               f"{max(Es)/(2*fl['VF2']*o['esc1']):.3f}")
    except FileNotFoundError:
        pass
    pr(f'[{time.time() - t0:.0f}s]')
    json.dump(rows, open(os.path.join(AF.CACHE, 'chi.json'), 'w'), default=float)
    open(os.path.join(AF.CACHE, 'chi.txt'), 'w').write('\n'.join(OUT) + '\n')


if __name__ == '__main__':
    main()
