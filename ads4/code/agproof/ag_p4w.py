"""agproof/ag_p4w.py (AG_PROOF.md §4): Lemma P4W's tolerances under WEAKER / RE-SHAPED decay premises than AG'.

P4_WINDOW's (AG'): |g_pm(p +- i y)| <= M(p) exp(-kappa_AG'(p, y)), kappa_AG' = 0.97 Lam beta/sqrt(1-x^2) for y <= 15/kappa, 0.97 Lam beta
beyond.  Here the edge (p = 1) and top majorants are replaced by one-parameter families (theta = the fraction of the rate assumed):
  'Fa'  : |g_pm(1 +- i y)| <= M(1) |flat_+(1 + i y)/flat_+(1)|^theta      (Lemma F flat wave, amplitude included; no threshold)
  'Fp'  : |g_pm(1 +- i y)| <= M(1) exp(-theta (J+1/2) Im[th(1+iy) - th(1)])  (Lemma F phase only)
  'U'   : |g_pm(1 +- i y)| <= M(1) exp(-theta Lam beta y)                     (uniform rate; no local enhancement at all)
  'Lt'  : P4_WINDOW's local rate with theta in place of 0.97 and NO threshold (rate theta Lam beta/sqrt(1-x1^2) on the whole side)
The top side uses the same family at height T = 25/(0.97 Lam beta) (min over p in [p_w, 1]); the boundary side at p_w and the window
terms are P4_WINDOW's (computed from the window series, not premised) and are left unchanged.  The curvature waves (P - J0)_pm keep
their explicit uniform rate (term A).
Everything else is p4win/p4w_states.analyse and the p4w_L loop (read-only imports; their contour_bound is replaced in this process only).
Output per state and family: the tolerated edge multiple Lam_AG (plain scalar envelope / +- transferred), and theta_min, the smallest
theta for which Lam_AG = 1 is tolerated.   -> cache/p4w.txt, cache/p4w.json
"""
import os, sys, json, time
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import mpmath as mp
import ag_flat as AF
import p4w_common as C          # read-only
import p4w_contour as K         # read-only
import p4w_states as PS         # read-only
import p4w_L as PL              # read-only

L = C.L
_orig_contour = K.contour_bound
FAM = {'name': 'Lt', 'theta': 0.97}
_prof_cache = {}
OUT = []
pr = lambda s: (print(s, flush=True), OUT.append(s))


def _yrule(kl, T):
    br = sorted(set([0.0] + [min(T, c/kl) for c in (1, 2, 4, 6, 9, 13, 20, 30, 45, 70)] + [T]))
    xs, ws = np.polynomial.legendre.leggauss(20)
    ys, wy = [], []
    for a, b in zip(br[:-1], br[1:]):
        if b > a:
            ys += list(a + (b - a)*(xs + 1)/2); wy += list(ws*(b - a)/2)
    return np.array(ys), np.array(wy)


def profiles(beta, mfl, T=None):
    """flat-wave profiles on the right side p = 1 and the top (exponents, theta = 1); rectangle height T (p-units)."""
    key = (round(beta, 12), round(mfl, 12), None if T is None else round(T, 12))
    if key in _prof_cache:
        return _prof_cache[key]
    m = mfl*L; Jh = m*np.sinh(beta/2)                # J + 1/2
    J = int(round(Jh - 0.5))
    x1 = 1/mfl; k0 = 0.97*beta*L; kl = k0/np.sqrt(1 - x1*x1)
    T = 25.0/k0 if T is None else T
    ys, wy = _yrule(kl, T)
    th = lambda z: complex(2*mp.asin(mp.mpc(z)*L/m))
    th1 = th(1.0)
    phase = np.array([Jh*(th(1 + 1j*y).imag - th1.imag) for y in ys])
    f1 = abs(AF.flatplus(J, m, L))
    amp = np.array([-float(mp.log(abs(AF.flatplus(J, m, L*(1 + 1j*y)))/f1)) for y in ys])
    pt = np.linspace(0, 1, 201)
    top_phase = np.array([Jh*th(p + 1j*T).imag for p in pt])            # Im th(p) = 0 on the real axis
    out = dict(ys=ys, wy=wy, phase=phase, amp=amp, T=T, kl=kl, k0=k0, pt=pt, top_phase=top_phase, J=J)
    _prof_cache[key] = out
    return out


def expo(pf, fam, theta):
    ys = pf['ys']
    if fam == 'Fa':
        return theta*pf['amp']
    if fam == 'Fp':
        return theta*pf['phase']
    if fam == 'U':
        return theta*beta_glob[0]*L*ys
    if fam == 'Lt':
        return theta/0.97*pf['kl']*ys
    raise ValueError(fam)


beta_glob = [0.0]


def contour_bound(beta, mfl, pw, Mw, Mfun, mode='uniform', k=2, u=1.0, Aunif=0.0):
    """Lemma P4W's (bd, edge, top) with the edge/top majorant of family FAM (mode 'local'); mode 'uniform' -> original."""
    if mode != 'local':
        return _orig_contour(beta, mfl, pw, Mw, Mfun, mode, k=k, u=u, Aunif=Aunif)
    bd, _, _ = _orig_contour(beta, mfl, pw, Mw, Mfun, 'local', k=k, u=u, Aunif=Aunif)     # boundary side: unchanged (computed)
    w = 2.0 if k == 2 else u
    fam, th = FAM['name'], FAM['theta']
    # rectangle height scaled so that the assumed decay at the top is still e^{-25 (theta/0.97)^0}: T = 25/(theta Lam beta)
    pf = profiles(beta, mfl, 25.0/(th*beta*L)); beta_glob[0] = beta
    ps = np.abs(K.psi_c(k, 1 + 1j*pf['ys']))
    Vp = float((ps*np.exp(-expo(pf, fam, th))*pf['wy']).sum())
    VA = K.V(k, 1.0, 0.97*beta*L, pf['T'])                     # curvature waves, explicit uniform rate
    M1 = float(Mfun(np.array([1.0]))[0])
    ed = 2*w*(Vp*(M1 + Aunif) + VA*Aunif)
    # top: |g_+(p + iT)| <= M_top(p) exp(-theta * Phi(p, T)); Phi from the flat phase (Fa, Fp, Lt) or theta beta T (U)
    x, wq = np.polynomial.legendre.leggauss(400); p = pw + (1 - pw)*(x + 1)/2; wq = wq*(1 - pw)/2
    if fam == 'U':
        et = np.exp(-th*beta*L*pf['T'])*np.ones_like(p)
    else:
        et = np.exp(-th*np.interp(p, pf['pt'], pf['top_phase']))
    tp = 2*w*float((np.abs(K.psi_c(k, p + 1j*pf['T']))*Mfun(p)*et*wq).sum())
    return bd, ed, tp


K.contour_bound = contour_bound        # this process only


def run_W(thetas, fams):
    runs = PS.PJ.load_runs(); ksh = PS.PJ.load_runs('ksh_*_*.json')
    states = json.load(open(os.path.join(C.CODE, 'p2con', 'cache', 'jstar_bh.json')))['states']
    states = sorted(states, key=lambda z: (z['J'] % 2, z['J'], z['mfl']))
    res = {}
    for fam in fams:
        for th in thetas:
            FAM['name'] = fam; FAM['theta'] = th
            rows = []
            for s in states:
                o = PS.analyse(s, runs, ksh)
                Lc = o['local']
                rows.append(dict(J=o['J'], mfl=o['mfl'], plain=Lc['LamAG_tol_plain'], tr=Lc['LamAG_tol'], fixed=Lc['fixed']/o['H'],
                                 d_plain=Lc['d1p']/o['H'], d_tr=Lc['d1']/o['H']))
            res[(fam, th)] = rows
            mp_ = min(r['plain'] for r in rows); mt = min(r['tr'] for r in rows)
            bad_p = [(r['J'], round(r['mfl'], 3)) for r in rows if r['plain'] < 1]
            bad_t = [(r['J'], round(r['mfl'], 3)) for r in rows if r['tr'] < 1]
            jp = min([r['J'] for r in rows if all(q['plain'] >= 1 for q in rows if q['J'] >= r['J'])] or [None])
            jt = min([r['J'] for r in rows if all(q['tr'] >= 1 for q in rows if q['J'] >= r['J'])] or [None])
            pr(f"  W47 family {fam:2s} theta = {th:.2f}: min Lam_AG tol plain {mp_:9.3g} / transferred {mt:9.3g}; sampled J* plain {jp}, "
               f"transferred {jt}; failing (plain) {bad_p[:6]}{'...' if len(bad_p) > 6 else ''}; failing (tr) {bad_t[:8]}{'...' if len(bad_t) > 8 else ''}")
    return res


def run_L(thetas, fams, La=3.0):
    res = {}
    for fam in fams:
        for th in thetas:
            FAM['name'] = fam; FAM['theta'] = th
            rows = []
            for lab, J, mtab, Fn, Fc in PL.ROWS:
                mfl = PL.mfl_exact(lab, J, mtab); m = mfl*L; beta = C.beta_of(J, mfl); u = 1/mfl**2
                nuw = 3.0/beta; pw = nuw/L
                nu = np.linspace(0, nuw, 4001)[1:]; dn = nu[1] - nu[0]
                P, dPb = C.Pd(beta, nu); P = P.real
                bb = (2*J + 1)/mfl*nu/L
                from scipy.special import j0
                B2 = (C.T(1, beta, nu, P, dPb)[0] - nu**2*C.gflat(1, bb))/m**2
                S4 = (C.T(2, beta, nu, P, dPb)[0] - nu**4*(C.gflat(2, bb) + (3/16)*j0(bb)))/m**4
                ps2 = C.psi2(nu/L); ps3 = C.psi3(nu/L)
                act2 = lambda f: 2*float((ps2*f).sum())*dn/L
                w2 = abs(act2(B2)) + u*float((np.abs(ps3)*np.abs(B2)).sum())*dn/L
                w4 = abs(act2(S4)) + u*float((np.abs(ps3)*np.abs(S4)).sum())*dn/L
                nn = np.linspace(0.7*nuw, 1.3*nuw, 41); Pn, dPn = C.Pd(beta, nn); bn = (2*J + 1)/mfl*nn/L
                B2n = (C.T(1, beta, nn, Pn.real, dPn)[0] - nn**2*C.gflat(1, bn))/m**2
                S4n = (C.T(2, beta, nn, Pn.real, dPn)[0] - nn**4*(C.gflat(2, bn) + (3/16)*j0(bn)))/m**4
                Mw = 2*float(np.max(np.abs(La*B2n + S4n)))
                env = K.e_sc(beta, mfl, hel=-2.0)
                Au = K.A_curv(beta, J, mfl)
                bd2, ed2, tp2 = K.contour_bound(beta, mfl, pw, Mw, env, 'local', k=2, Aunif=Au)
                bd3, ed3, tp3 = K.contour_bound(beta, mfl, pw, Mw, env, 'local', k=3, u=u, Aunif=Au)
                d = La*w2 + w4 + bd2 + bd3 + ed2 + tp2 + ed3 + tp3
                edge = ed2 + ed3 + tp2 + tp3
                tol = (Fc/2 - (d - edge))/edge            # tolerated multiple of the transferred scalar edge envelope
                rows.append(dict(label=lab, J=J, mfl=mfl, d=d/(Fc/2), edge=edge/(Fc/2), tol=tol))
            res[(fam, th)] = rows
            pr(f"  L' family {fam:2s} theta = {th:.2f}: max d/(F_cert/2) = {max(r['d'] for r in rows):.3f}; max edge+top share "
               f"{max(r['edge'] for r in rows):.1e}; min edge tolerance (multiple of the +- transferred envelope) {min(r['tol'] for r in rows):.3g} "
               f"(at J = {min(rows, key=lambda r: r['tol'])['J']}, {min(rows, key=lambda r: r['tol'])['label']})")
    return res


def main():
    t0 = time.time()
    pr(__doc__.split('->')[0])
    fams = ['Lt', 'Fa', 'U']
    thetas = [0.97, 0.8, 0.6, 0.45]
    pr('\n== Theorem W, 47 sampled states (c = 3, theta_H = 1); P4_WINDOW reference (AG\', T1 = 15): min tol plain 114, transferred 0.287')
    rw = run_W(thetas, fams)
    pr('\n== Theorem L\' states (Lam_amp = 3, +- transferred envelope, certified half-floor); P4_WINDOW reference: max d = 0.086')
    rl = run_L([0.97, 0.6, 0.4], fams)
    pr(f'[{time.time() - t0:.0f}s]')
    json.dump(dict(W={f'{k[0]}|{k[1]}': v for k, v in rw.items()}, L={f'{k[0]}|{k[1]}': v for k, v in rl.items()}),
              open(os.path.join(AF.CACHE, 'p4w.json'), 'w'), default=float)
    open(os.path.join(AF.CACHE, 'p4w.txt'), 'w').write('\n'.join(OUT) + '\n')


if __name__ == '__main__':
    main()
