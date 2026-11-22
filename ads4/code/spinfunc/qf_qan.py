"""
qf_qan (RESULT_Q_FWD_CLOSE.md): Result Q chain (logic copied from bh_qan4 / po_qan5, REFCORR = 1) for the channels measured here.
  Gamma_exact = m^2 (c T_ref I_st/I_ref - T_st) + gamma_h (m/m_ref)^2 sum_j R_j (nu - ij)^2,
  Gamma_K     = Gamma_exact - gamma_h sum_j R_j (nu - ij)^2 - const N^_4  (Gamma_K(0.02) = 0),
  s_d, s_v = derivative / value projections per gauge (bf_qan.weights: H0 H2 H3 H4 SB XP), truncation = |s(K) - s(K-1)|,
  FF/N at the state from the Stage-6 HYB model (bg_qan.ffn), FF -> FF - a 2u^3 s_d - c (2u^3/Lam^2) s_v.
Modes
  oo    D J K Dg : parity-odd structure at even J (OO channel), 1 x 1 (qf_eeq odd cache); gamma = gamma_++ = 0.75, HYB row '++'.
                   Also the ++ image at the SAME Delta (bg_eeq cache, or qf_eeq even cache) through the SAME chain, and the
                   OO - ++ comparison: per harmonic, s_d, and the reference ratio I_OO,ref/I_++,ref (order-0 lemma: nu-constant).
  pmodd D J K Dg : parity-odd +- at odd J (as po_qan5), gamma = gamma_+- = 20.75.
  even  D J K Dg : parity-even 2 x 2 at J >= 4 in the order-0 helicity basis T0(D) (qf_helbasis), ++ and +- (as bh_qan4) + mixing.
Usage: python3 qf_qan.py mode D J K Dg  -> appends cache/qf_qan.txt, pickle cache/qf_qan_{mode}_D{D}_J{J}_Dg{Dg}.pkl
"""
import os, sys, pickle
os.environ.setdefault('L1_DEG', '6')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
RQ = os.path.join(HERE, '..', 'resultq'); sys.path.insert(0, RQ)
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import bf_qan as BQ
import bg_eeq as EQ
import bg_qan as GQ
import qf_eeq as QE

GAM = {'++': 0.75, '+-': 20.75, 'OO': 0.75}
HYBROW = {'++': '++', '+-': '+-', 'OO': '++'}
XSMIN = float(os.environ.get('QF_XSMIN', '0.15'))     # helicity-basis x range (l8 rule: 0.15); 0.25 = systematic check
N = 1024
NUP = (0.02, 0.5, 1.02, 1.5, 2.02, 2.5, 3.02, 4.02, 5.02, 6.02, 8.02, 10.02, 13.02, 16.02, 20.02, 30.02, 40.02)


def tcache(nus, eta, R, key):
    for pre in ('bh_qan_T_', 'po_qan_T_', 'bg_qan_T_', 'qf_qan_T_'):
        fn = os.path.join(HERE, 'cache', pre + key + '.pkl')
        if os.path.exists(fn):
            d = pickle.load(open(fn, 'rb'))
            if len(d) == len(nus):
                return d
    v = BQ.T(nus, eta, R)
    pickle.dump(v, open(os.path.join(HERE, 'cache', 'qf_qan_T_' + key + '.pkl'), 'wb'))
    return v


def base(D, J):
    nus = EQ.NUS[:N]
    m, eta = BQ.mass(D, J); mref, etaref = BQ.mass(EQ.DREF, J)
    R = BQ.Rj(nus)
    Tst = tcache(nus, eta, R, 'D%g_J%d' % (D, J)); Tref = tcache(nus, etaref, R, 'D%g_J%d' % (EQ.DREF, J))
    NH = sum(R.values()); cov = sum(R[j]*(nus - 1j*j)**2 for j in R)
    AR = sum(np.abs(v) for v in R.values())
    return dict(nus=nus, m=m, mref=mref, Tst=Tst, Tref=Tref, NH=NH, cov=cov, AR=AR)


def chain(B, Is, Ir, h):
    tot = B['Tref']*Is/Ir
    cc = B['Tst'][0]/tot[0]
    m = B['m']
    Gex = m*m*(cc*tot - B['Tst']) + GAM[h]*(m/B['mref'])**2*B['cov']
    GK = Gex - GAM[h]*B['cov']
    GK = GK - (GK[0]/B['NH'][0])*B['NH']
    return dict(GK=GK, Gex=Gex, ratio=Is/Ir, c=cc)


def images_1x1(st, rf, Kt):
    return EQ.value(st, 'EE4', 0, 0, Kt)[0][:N], EQ.value(rf, 'EE4', 0, 0)[0][:N]


def helmats(rec, D, J, Kt=None):

    import qf_helbasis as HB
    binv = shadow_inv(D, J)
    Ti = HB.basis(J, D, xs_min=XSMIN)['Ti']
    I = np.array([[EQ.value(rec, 'EE4', a, c, Kt)[0][:N] for c in range(2)] for a in range(2)])
    M = np.einsum('acn,bc->nab', I, binv)
    return np.einsum('ia,nab,jb->nij', Ti, M, Ti)


def shadow_inv(D, J):
    fn = os.path.join(HERE, 'cache', 'qf_shadow.pkl')
    db = pickle.load(open(fn, 'rb')) if os.path.exists(fn) else {}
    key = (round(float(D), 6), J)
    if key not in db:
        bh = os.path.join(HERE, 'cache', 'bh_shadow.pkl')
        bdb = pickle.load(open(bh, 'rb')) if os.path.exists(bh) else {}
        if key in bdb:
            db[key] = bdb[key]
        else:
            import l6_spin_action as L6
            b, logK, resid = L6.shadow(float(D), J)
            db[key] = dict(b=b, logK=logK, resid=resid)
        pickle.dump(db, open(fn, 'wb'))
    return np.linalg.inv(db[key]['b'])


def report(pr, B, rs, K, h, J, u, Dg, W, tag=''):
    r = rs[K]
    idx = [int(np.argmin(np.abs(B['nus'] - v))) for v in NUP]
    AR = B['AR']
    pr(f'  [{tag or h}]  nu: Gamma_K/sum|R_j| (K) [a^K term]  |  Gamma_exact/sum|R_j|  |  Re Gex/Re cov')
    for i in idx:
        g = r[h]['GK'][i]; gp = rs[K - 1][h]['GK'][i]
        pr(f'     {B["nus"][i]:5.2f}: {g.real/AR[i]:+10.3f} [{(g - gp).real/AR[i]:+8.3f}]  | {r[h]["Gex"][i].real/AR[i]:+10.3f}'
           f'  | {r[h]["Gex"][i].real/B["cov"][i].real if abs(B["cov"][i].real) > 0 else 0:+8.3f}')
    P = {Kt: GQ.project(dict(nus=B['nus'], GK=rs[Kt][h]['GK']), W) for Kt in range(K + 1)}
    P18 = GQ.project(dict(nus=B['nus'], GK=r[h]['GK']), W, 18.0)
    out = dict(P=P)
    pr(f'  [{tag or h}] projections (nu < 41; [nu < 18]): s_d +- trunc | s_v +- trunc | s_d by order k = 0..K | FF/N (HYB row {HYBROW[h]}): K=0 -> measured')
    for g in W:
        sd, sv = P[K][g]; ed = abs(sd - P[K - 1][g][0]); ev = abs(sv - P[K - 1][g][1])
        f = GQ.ffn(HYBROW[h], J, u, Dg, sd, sv)
        pr(f'     {g:3s}: s_d = {sd:+9.4f} +- {ed:.4f} [{P18[g][0]:+8.4f}] | s_v = {sv:+9.3f} +- {ev:.3f} | '
           + ' '.join(f'{P[k][g][0] - (P[k-1][g][0] if k else 0):+.4f}' for k in range(K + 1)) + ' | '
           + '  '.join(f'd={dl}: {f[dl][0]:+.4f} -> {f[dl][1]:+.4f} (+- {f[dl][2]*ed:.4f})' for dl in f))
        out[g] = dict(sd=sd, ed=ed, sv=sv, ev=ev, ff=f)
    w = 1/AR; sel = B['nus'] < 12.0
    A = np.array([B['cov'].real[sel]*w[sel], B['NH'].real[sel]*w[sel]]).T
    co, *_ = np.linalg.lstsq(A, r[h]['GK'].real[sel]*w[sel], rcond=None)
    rem = r[h]['GK'] - co[0]*B['cov'] - co[1]*B['NH']
    Pr = GQ.project(dict(nus=B['nus'], GK=rem), W)
    pr(f'  [{tag or h}] covariant part dgamma = {co[0]:+.4f} (= {co[0]/GAM[h]*100:+.2f} % of gamma_h); remainder s_d: '
       + ' '.join(f'{g}:{Pr[g][0]:+.4f}' for g in W))
    out['cov'] = dict(dgamma=co[0], rem={g: Pr[g][0] for g in W})
    return out


if __name__ == '__main__':
    mode, D, J, K, Dg = sys.argv[1], float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    W = BQ.weights()
    lines = []

    def pr(s=''):
        print(s, flush=True); lines.append(s)
    B = base(D, J)
    Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0)); u = (Lam/B['m'])**2
    res = {}; rs = {}
    pr(f'# qf_qan {mode} (helicity-basis xs_min {XSMIN}): EE k=4, Delta={D:.4f} J={J} m={B["m"]:.3f} = {B["m"]/Lam:.4f} Lam(Dg={Dg}), K={K}, REFCORR=1')
    if mode in ('oo', 'pmodd'):
        h = 'OO' if mode == 'oo' else '+-'
        st = QE.load('odd', D, J, K); rf = QE.load('odd', EQ.DREF, J, 2)
        for Kt in range(K + 1):
            rs[Kt] = {h: chain(B, *images_1x1(st, rf, Kt), h)}
        res[h] = report(pr, B, rs, K, h, J, u, Dg, W)
        if mode == 'oo':
            # ++ at the same Delta through the same chain
            pp = None
            p1 = EQ.fn(D, J, K)
            if os.path.exists(p1):
                db = pickle.load(open(p1, 'rb')); pst, prf = db['st'], db['ref']
                pp = 'bg_eeq'
            else:
                try:
                    pst = QE.load('even', D, J, K); prf = QE.load('even', EQ.DREF, J, 2); pp = 'qf_eeq even'
                except FileNotFoundError:
                    pst = None
            if pst is not None:
                for Kt in range(K + 1):
                    rs[Kt]['++'] = chain(B, *images_1x1(pst, prf, Kt), '++')
                res['++'] = report(pr, B, rs, K, '++', J, u, Dg, W, tag=f'++ same Delta ({pp})')
                # comparison
                idx = [int(np.argmin(np.abs(B['nus'] - v))) for v in NUP]
                AR = B['AR']
                Ior = EQ.value(rf, 'EE4', 0, 0)[0][:N]; Ipr = EQ.value(prf, 'EE4', 0, 0)[0][:N]
                rr = Ior/Ipr; rr = rr/rr[0]
                pr('  OO vs ++ (helicity lemma test, k = 4):')
                pr('     reference ratio I_OO,ref/I_++,ref / (value at nu=0.02) [order-0 lemma: = 1]: '
                   + ' '.join(f'{rr[i].real:+.7f}{rr[i].imag:+.1e}i' for i in idx[::2]))
                pr('     max |ratio - 1| (nu < 41): %.2e;  m_ref^2 x that: %.3f c2' % (np.max(np.abs(rr - 1)), B['mref']**2*np.max(np.abs(rr - 1))))
                for Kt in (K,):
                    d = rs[Kt]['OO']['GK'] - rs[Kt]['++']['GK']
                    pr(f'     (Gamma_K^OO - Gamma_K^++)/sum|R_j| (K={Kt}): ' + ' '.join(f'{(d[i]/AR[i]).real:+.1e}' for i in idx))
                    pr(f'     |Gamma_K^OO - Gamma_K^++| max on nu<2.5 / 2.5-18 / 18-41: '
                       + ' / '.join(f'{np.max(np.abs(d[(B["nus"] >= lo) & (B["nus"] < hi)])):.3f}' for lo, hi in ((0, 2.5), (2.5, 18), (18, 41)))
                       + ' c2;  per |R|: ' + ' / '.join(f'{np.max(np.abs(d/AR)[(B["nus"] >= lo) & (B["nus"] < hi)]):.1e}'
                                                         for lo, hi in ((0, 2.5), (2.5, 18), (18, 41))))
                dd = {Kt: rs[Kt]['OO']['GK'] - rs[Kt]['++']['GK'] for Kt in range(K + 1)}
                Pd = {Kt: GQ.project(dict(nus=B['nus'], GK=dd[Kt]), W) for Kt in range(K + 1)}
                pr('     s_d(OO) - s_d(++) per gauge [by order k = 0..K]: ' + '; '.join(
                    f'{g}: {Pd[K][g][0]:+.1e} [' + ' '.join(f'{Pd[k][g][0] - (Pd[k-1][g][0] if k else 0):+.1e}' for k in range(K + 1)) + ']' for g in W))
                res['diff'] = {g: Pd[K][g][0] for g in W}
                res['diff_orders'] = {g: [Pd[k][g][0] - (Pd[k-1][g][0] if k else 0) for k in range(K + 1)] for g in W}
                res['refratio_dev'] = float(np.max(np.abs(rr - 1)))
    elif mode == 'even':
        import qf_helbasis as HB
        st = QE.load('even', D, J, K); rf = QE.load('even', EQ.DREF, J, 2)
        Bh = HB.basis(J, D, xs_min=XSMIN); Br = HB.basis(J, EQ.DREF, xs_min=XSMIN)
        pr(f'  helicity basis T0(D): sv3/sv1 = {Bh["sv3"]:.1e}, order-0 offdiag {Bh["offdiag0"]:.1e}, flatness {Bh["flat"][0]:.2e} / {Bh["flat"][1]:.2e};'
           f' reference: {Br["sv3"]:.1e}, {Br["offdiag0"]:.1e}')
        Mr = helmats(rf, EQ.DREF, J)
        for Kt in range(K + 1):
            Ms = helmats(st, D, J, Kt)
            rs[Kt] = {}
            cs = {}
            for i, h in enumerate(('++', '+-')):
                rs[Kt][h] = chain(B, Ms[:, i, i], Mr[:, i, i], h); cs[h] = rs[Kt][h]['c']
            m = B['m']
            off = m*m*np.sqrt(cs['++']*cs['+-'])*B['Tref']*0.5*(Ms[:, 0, 1] + Ms[:, 1, 0])/np.sqrt(Mr[:, 0, 0]*Mr[:, 1, 1])
            rs[Kt]['off'] = dict(GK=off, kappa=0.5*(Ms[:, 0, 1] + Ms[:, 1, 0])/np.sqrt(Ms[:, 0, 0]*Ms[:, 1, 1]))
            if Kt == K:
                rref = Mr[:, 1, 1]/Mr[:, 0, 0]
                pr(f'  reference: |offdiag|/sqrt(diag) max {np.max(np.abs(0.5*(Mr[:,0,1]+Mr[:,1,0])/np.sqrt(Mr[:,0,0]*Mr[:,1,1]))):.1e}; '
                   f'I_+-/I_++ at nu = 0.02/10/40: {rref[0].real:.6f} {rref[250].real:.6f} {rref[1000].real:.6f}')
        for h in ('++', '+-'):
            res[h] = report(pr, B, rs, K, h, J, u, Dg, W)
        idx = [int(np.argmin(np.abs(B['nus'] - v))) for v in NUP]
        Po = {Kt: GQ.project(dict(nus=B['nus'], GK=rs[Kt]['off']['GK']), W) for Kt in (K - 1, K)}
        pr('  mixing: m^2 kappa: ' + ' '.join(f'{(B["m"]**2*rs[K]["off"]["kappa"][i]).real:+.2f}' for i in idx[:12])
           + ' | s_off per gauge: ' + ' '.join(f'{g}:{Po[K][g][0]:+.4f}(+-{abs(Po[K][g][0] - Po[K-1][g][0]):.4f})' for g in W))
        res['off'] = {g: Po[K][g] for g in W}
        d = rs[K]['++']['GK'] + rs[K]['+-']['GK']
        Pe = GQ.project(dict(nus=B['nus'], GK=0.5*d), W)
        pr('  helicity-even part (Gamma_++ + Gamma_+-)/2: s_d ' + ' '.join(f'{g}:{Pe[g][0]:+.4f}' for g in W)
           + ' | per |R|: ' + ' '.join(f'{(0.5*d[i]/B["AR"][i]).real:+.2f}' for i in idx[:12]))
    open(os.path.join(HERE, 'cache', 'qf_qan.txt'), 'a').write('\n'.join(lines) + '\n')
    pickle.dump(dict(rs=rs, res=res, D=D, J=J, Dg=Dg, m=B['m'], u=u, mode=mode, K=K),
                open(os.path.join(HERE, 'cache', 'qf_qan_%s_D%g_J%d_Dg%d%s.pkl' % (mode, D, J, Dg, '' if XSMIN == 0.15 else '_xs%g' % XSMIN)), 'wb'))
