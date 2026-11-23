"""p4win/p4w_states.py (P4_WINDOW.md §4): Lemma P4W / Theorem W-win on Theorem W's 47 sample states.

Per state and entry e of M_Omega (even J: (+-,+-), (++,++), (++,+-); odd J: (+-,+-)), with G_e := exact - hybrid (de-shifted, per unit
leading image) and the Regge window nu < nu_w = c/beta (c = 3 unless stated):
  window, computed  : w2 = |2 sum_{nu<nu_w} psi2 f2bh dnu/Lam| / m^2      f2bh = F2/c (block + kernel O(r), p2con caches) minus the hybrid's
                                                                           own flat O(y) tower nu^2 c1_h(bhat p) (pc_jstar.flat_c1; 0 off-diag.)
                      w3 = u sum_{nu<nu_w} |psi3| |f2bh| dnu/Lam / m^2    (k = 3, same structures, sign-free: INDEP_CHECK2 convention)
                      w4 = |window cost of S4/m^4| (scalar seed, beyond-hybrid O(1/m^4) = T_2 - nu^4 g_2), + S6/m^6, S8/m^8 (ratio test)
  window, premise   : K4w * c4w,  c4w = window cost of |eta^2 nu^4 P| + |eta sinh(beta) nu^4 dP/dbeta|, per 1/m^4 (spinning O(1/m^4),
                      kernel O(r^2): not computed; K4w is the tolerated coefficient)
  monopole          : 2|D||G(i/2)|: O(1/m^2) part 0 (computed, P2_CONTRACT §1.7); O(1/m^4) tolerance Q4_tol = theta H m^4/(2|D|)
  complement (AG')  : contour bound of Lemma P4W (p4w_contour): bd = boundary term at p_w with M_w = 2 sup_{[0.7,1.3] nu_w} |f2bh/m^2 + S4/m^4|
                      (computed), and Lam_AG x (edge + top) with the scalar-calibrated envelope e_sc (p4w_check: valid on the exact scalar
                      with margin 1.5), local rate kappa = 0.97 Lam beta / sqrt(1 - x^2) ('local') or uniform 0.97 Lam beta ('uniform').
Criterion (Theorem W, sign-free 2x2): d := max_e [w2 + w3 + w4 + bd2 + bd3] + Lam_AG (edge + top)_{k=2,3} <= theta H.
Output: tolerated Lam_AG (local / uniform), tolerated K4w, and the split of d at Lam_AG = 1.  -> cache/states.txt, cache/states.json
"""
import os, sys, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.special import j0
import p4w_common as C, p4w_contour as K
sys.path.insert(0, os.path.join(C.CODE, 'p2con'))
import pc_jstar as PJ          # read-only (load_runs, flat_c1; through it pc_run.Feval, cfac)
Rn = PJ.Rn

OUT = []
pr = lambda s: (print(s, flush=True), OUT.append(s))


def f2bh(e_xs, beta, J, mfl, h, nu):
    f = 0
    for x in e_xs:
        f = f + (Rn.Feval(np.array(x), beta, nu)[0]/Rn.cfac(nu)).real
    if h is not None:
        f = f - nu**2*PJ.flat_c1(h, (2*J + 1)/mfl*nu/C.L)
    return f


def scal4(beta, J, mfl, m, nu):
    P, dPb = C.Pd(beta, nu)
    bb = (2*J + 1)/mfl*nu/C.L
    S = {k: (C.T(k, beta, nu, P.real, dPb)[0] - nu**(2*k)*C.gflat(k, bb))/m**(2*k) for k in (2, 3, 4)}
    S[2] = S[2] - (3/16)*nu**4*j0(bb)/m**4        # the scalar hybrid's flat O(y^2) term is g_2 + (3/16) J0 (p4w_316)
    return S, P.real, dPb.real


def analyse(s, runs, ksh, c=3.0, theta=1.0):
    J, mfl, m, beta, H = s['J'], s['mfl'], s['m'], s['beta'], s['H']
    u = 1/mfl**2
    key = (J, round(mfl, 6)); r = runs[key]; kk = ksh[key]
    e = r['ent']; ke = kk['ent']
    if J % 2 == 0:
        pairs = {('+-', '+-'): ('same|same', '+-'), ('++', '++'): ('opp|opp', '++'), ('++', '+-'): ('same|opp', None)}
    else:
        pairs = {('+-', '+-'): ('odd|odd', '+-')}
    nuw = c/beta; pw = nuw/C.L
    sel = C.W.NU < nuw
    nu = C.W.NU[sel]; ps2 = C.W.PS[sel]; ps3 = C.psi3(nu/C.L); dn = C.W.DN
    S, P, dPb = scal4(beta, J, mfl, m, nu)
    act2 = lambda f: 2*float((ps2*f).sum())*dn/C.L
    w4 = abs(act2(S[2])); w6 = abs(act2(S[3])); w8 = abs(act2(S[4]))
    q = max(w6/max(w4, 1e-300), w8/max(w6, 1e-300))
    tail = w8*q/(1 - q) if q < 1 else np.inf
    Iw = (2*float(np.abs(ps2).sum()) + u*float(np.abs(ps3).sum()))*dn/C.L      # 2 int_0^{p_w} |psi2| + u int |psi3|
    supS4 = float(np.max(np.abs(S[2])))
    c4w = (abs(act2(np.cosh(beta)**2*nu**4*P)) + abs(act2(np.cosh(beta)*np.sinh(beta)*nu**4*dPb)))/m**4
    nn = np.linspace(0.7*nuw, 1.3*nuw, 41)
    Sn, Pn, dPn = scal4(beta, J, mfl, m, nn)
    m4w = 2*float(np.max(np.cosh(beta)**2*nn**4*np.abs(Pn) + np.cosh(beta)*np.sinh(beta)*nn**4*np.abs(dPn)))/m**4
    ent = {}
    for k_, (v, h) in pairs.items():
        xs = [e[v + '|full']['x'], ke[v]['x']]
        f = f2bh(xs, beta, J, mfl, h, nu)
        w2 = abs(act2(f))/m**2
        w3 = u*float((np.abs(ps3)*np.abs(f)).sum())*dn/C.L/m**2
        fn = f2bh(xs, beta, J, mfl, h, nn)
        Mw = 2*float(np.max(np.abs(fn/m**2 + Sn[2])))
        ent['%s,%s' % k_] = dict(w2=w2, w3=w3, Mw=Mw)
    Mw = max(v['Mw'] for v in ent.values())
    out = dict(J=J, mfl=mfl, m=m, beta=beta, H=H, c=c, nuw=nuw, pw=pw, ent=ent, w4=w4, w6=w6, w8=w8, q=q, tail=tail, c4w=c4w, Mw=Mw,
               x1=1/mfl)
    HEL = {'+-,+-': -2.0, '++,++': 2.0, '++,+-': 0.0}
    x1 = 1/mfl; E2h = (2 - x1*x1)/(1 - x1*x1)**2/2
    Au = {ek: (0.0 if ek == '++,+-' else K.A_curv(beta, J, mfl)) for ek in ent}
    for mode in ('local', 'uniform'):
        win = max(v['w2'] + v['w3'] for v in ent.values()) + w4 + w6 + w8 + (tail if np.isfinite(tail) else 0)
        per = {}
        for ek in ent:
            zero = lambda p: 0.0*np.asarray(p, float)
            b0, e0, t0_ = K.contour_bound(beta, mfl, pw, Mw, zero, mode, k=2, Aunif=Au[ek])
            b0_, e0_, _ = K.contour_bound(beta, mfl, pw, Mw, zero, mode, k=3, u=u, Aunif=Au[ek])
            A_only = e0 + e0_                                  # edge cost of the known curvature waves alone (envelope M = 0)
            for var, hl in (('tr', HEL[ek]), ('pl', 0.0)):
                env = K.e_sc(beta, mfl, hel=hl)
                bd2, ed2, tp2 = K.contour_bound(beta, mfl, pw, Mw, env, mode, k=2, Aunif=Au[ek])
                bd3, ed3, tp3 = K.contour_bound(beta, mfl, pw, Mw, env, mode, k=3, u=u, Aunif=Au[ek])
                M1 = float(env(np.array([1.0]))[0])
                per[ek + '|' + var] = dict(bd=bd2 + bd3, A_only=A_only, ag_env=ed2 + ed3 - A_only + tp2 + tp3, ag=ed2 + tp2 + ed3 + tp3,
                                           slope=(ed2 + ed3 - A_only)/M1, M1=M1, env_hyb=E2h*(1 - x1*x1)**HEL[ek])
        bd = max(v['bd'] for v in per.values()); fixed = win + bd
        bdpm = bd/Mw
        trs = [v for k_, v in per.items() if k_.endswith('|tr')]; pls = [v for k_, v in per.items() if k_.endswith('|pl')]
        ag = max(v['ag'] for v in trs); agp = max(v['ag'] for v in pls)
        LamT = min((theta*H - fixed - v['A_only'])/v['ag_env'] for v in trs)
        LamP = min((theta*H - fixed - v['A_only'])/v['ag_env'] for v in pls)
        Medge = min((theta*H - fixed - v['A_only'])/v['slope'] for v in pls)        # tolerated edge amplitude M(1), any entry
        rho_tol = min((theta*H - fixed - v['A_only'])/v['slope']/v['env_hyb'] for v in trs)
        out[mode] = dict(bd=bd, win=win, fixed=fixed, ag=ag, agp=agp, per=per, LamAG_tol=LamT, LamAG_tol_plain=LamP,
                         A_only=max(v['A_only'] for v in pls),
                         K4w_tol=(theta*H - fixed)/(c4w + bdpm*m4w),
                         K4w_tol_AG=(theta*H - fixed - agp)/(c4w + bdpm*m4w), d1=fixed + ag, d1p=fixed + agp,
                         Medge_tol_c2=Medge*m**2, Msc_edge_c2=per['+-,+-|pl']['M1']*m**2, rho_tol_pm=rho_tol)
    out['Q4_tol'] = theta*H*m**4/(2*abs(C.D1000))
    out['Iw'] = Iw; out['supS4'] = supS4
    out['epsw_tol'] = (theta*H - out['local']['bd'] - out['local']['agp'] - max(v['w2'] + v['w3'] for v in ent.values()))/Iw
    return out


def main():
    t0 = time.time()
    pr(__doc__.split('Output')[0])
    runs = PJ.load_runs(); ksh = PJ.load_runs('ksh_*_*.json')
    states = json.load(open(os.path.join(C.CODE, 'p2con', 'cache', 'jstar_bh.json')))['states']
    p4b = {(q['J'], round(q['mfl'], 6)): q for q in json.load(open(os.path.join(C.CODE, 'p2con', 'cache', 'p4b.json')))}
    res = []
    pr('  J   m/L   beta  nu_w  |  H        | window (max entry: w2+w3 | O(1/m^4..8) scalar, ratio q) | bd (k=2+3) | AG edge+top at Lam_AG=1: '
       'plain / helicity-transferred | d(Lam_AG=1)/H plain / transferred | Lam_AG tol (local): plain / transferred (uniform rate: plain) | '
       'edge amplitude tol (c2) vs scalar m^2 e_sc(x1) | rho_tol(+-) | K4w tol (AG=0 / AG plain=1) | P2_CONTRACT full-range K4(R4) vs tol')
    for s in sorted(states, key=lambda z: (z['J'] % 2, z['J'], z['mfl'])):
        o = analyse(s, runs, ksh)
        L_ = o['local']; U = o['uniform']
        wmax = max(v['w2'] + v['w3'] for v in o['ent'].values())
        q = p4b.get((s['J'], round(s['mfl'], 6)), {})
        pr(f"  {o['J']:3d} {o['mfl']:.3f} {o['beta']:.4f} {o['nuw']:5.1f} | {o['H']:.2e} | {wmax/o['H']:.1e} | {(o['w4'] + o['w6'] + o['w8'])/o['H']:.1e} "
           f"(q {o['q']:.0e}) | {L_['bd']/o['H']:.2f} | {L_['agp']/o['H']:.1e} / {L_['ag']/o['H']:.1e} | {L_['d1p']/o['H']:.2f} / {L_['d1']/o['H']:.2f} | "
           f"{L_['LamAG_tol_plain']:7.3g} / {L_['LamAG_tol']:7.3g} ({U['LamAG_tol_plain']:.2g}) | {L_['Medge_tol_c2']:.1e} vs {L_['Msc_edge_c2']:.1e} | "
           f"{L_['rho_tol_pm']:.3g} | {L_['K4w_tol']:7.3g} / {L_['K4w_tol_AG']:7.3g} | {q.get('K_rem', np.nan):7.3g} vs {s.get('K4_tol', np.nan):6.3g}")
        o['K4_series_fullrange'] = q.get('K_rem'); o['K4_tol_W'] = s.get('K4_tol')
        res.append(o)
    pr('\nsummary (c = 3, theta = 1, local rate): min Lam_AG tol plain = %.3g, transferred = %.3g; min K4w tol (AG = 0) = %.3g; '
       'max (window + bd)/H = %.3g; max d(Lam_AG = 1, plain)/H = %.3g' % (
           min(o['local']['LamAG_tol_plain'] for o in res), min(o['local']['LamAG_tol'] for o in res), min(o['local']['K4w_tol'] for o in res),
           max(o['local']['fixed']/o['H'] for o in res), max(o['local']['d1p']/o['H'] for o in res)))
    pr('\nwindow sup-norm tolerance eps_w (|G - F2bh/m^2| <= eps_w on nu < nu_w, real axis only) vs sup |S4/m^4| (scalar O(1/m^4)):')
    pr('  ' + '; '.join(f"({o['J']}, {o['mfl']:.3f}): {o['epsw_tol']:.1e} vs {o['supS4']:.1e}" for o in res))
    pr('  min eps_w tol / sup|S4| = %.3g' % min(o['epsw_tol']/o['supS4'] for o in res))
    for thr in (1.0, 3.0, 10.0, 30.0):
        bad = sorted(set((o['J'], round(o['mfl'], 3)) for o in res if o['local']['LamAG_tol'] < thr))
        pr(f'  transferred envelope, Lam_AG = {thr:g}: failing sampled states (J, m/L): {bad}')
    # sensitivity to c
    sens = {}
    for c in (2.0, 5.0):
        lt = []
        for s in states:
            o = analyse(s, runs, ksh, c=c)
            lt.append((o['local']['LamAG_tol'], o['J'], o['local']['d1']/o['H'], o['mfl'], o['local']['fixed']/o['H'], o['local']['LamAG_tol_plain']))
        sens[c] = lt
        pr(f'  c = {c:g}: min Lam_AG tol (local, transferred) = {min(lt)[0]:.3g} at J = {min(lt)[1]}; plain {min(v[5] for v in lt):.3g}; max d(Lam_AG=1, transferred)/H = {max(v[2] for v in lt):.3g}; '
           f'max (window + bd)/H = {max(v[4] for v in lt):.3g}')
        pr('     per state (J, m/L): fixed/H, Lam_AG tol: ' + '; '.join(f'({v[1]}, {v[3]:.3f}): {v[4]:.2f}, {v[0]:.3g}' for v in sorted(lt, key=lambda v: (v[1] % 2, v[1], v[3]))))
    pr(f'[{time.time() - t0:.0f}s]')
    json.dump(dict(c3=res, sens={str(k): v for k, v in sens.items()}), open(os.path.join(C.CACHE, 'states%s.json' % os.environ.get('P4W_T1', '')), 'w'), default=float)
    open(os.path.join(C.CACHE, 'states%s.txt' % os.environ.get('P4W_T1', '')), 'w').write('\n'.join(OUT) + '\n')


if __name__ == '__main__':
    main()
