"""
bh_qan4 (BELOW_FOLD3.md §1): Result Q chain (bg_qan, copied logic) for a J >= 4 state with two conserved <TTO> structures, per helicity.
  M_{ab}(nu) = sum_c I_{ac}(nu) (b^-1)_{bc}       (bg prefactor matrix x inverse shadow matrix, as l6_spin_action.action: M = Ser binv^T)
  Mt(nu)     = Ti M(nu) Ti^T,  Ti from bh_helbasis (the project's order-0 Laplace helicity basis at the SAME Delta; reference: its own)
  I_h = Mt_hh;  Gamma_exact,h = m^2 (c_h T_ref I_h,st/I_h,ref - T_st);  Gamma_K,h = Gamma_exact,h - Gamma_HYB,h - const N^_4,
  Gamma_HYB,h = gamma_h (1 - [m/m_ref]^2 REFCORR) sum_j R_j (nu - ij)^2 ... see note:  the reference's own finite-J term gamma_h (m/m_ref)^2
  sum_j R_j (nu-ij)^2 is ADDED back (REFCORR = 1, default; = 0 reproduces the bg_qan convention, which neglects it: 2e-4 x gamma_h).
  s_d, s_v per gauge (bf_qan.weights), truncation = |s(K) - s(K-1)|; FF/N at the state (bg_qan.ffn, Stage-6 HYB).
  Off-diagonal (helicity mixing): Gamma_off = m^2 sqrt(c_++ c_+-) T_ref Mt_01,st / sqrt(I_++,ref I_+-,ref) and its derivative projection.
Usage: python3 bh_qan4.py Delta J K Dg [h-list]  (needs cache/bh_eeq_D{Delta}_J{J}_K{K}.pkl and the reference bh_eeq_D20000_J{J}_K2.pkl)
       -> appends cache/bh_qan4.txt, pickle cache/bh_qan4_D{Delta}_J{J}.pkl
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
import bh_eeq4 as B4
import bh_helbasis as HB

GAM = {'++': 0.75, '+-': 20.75}
REFCORR = float(os.environ.get('BH_REFCORR', '1'))
N = 1024


def tcache(nus, eta, R, key):
    fn = os.path.join(HERE, 'cache', 'bh_qan_T_%s.pkl' % key)
    if os.path.exists(fn):
        d = pickle.load(open(fn, 'rb'))
        if len(d) == len(nus):
            return d
    v = BQ.T(nus, eta, R)
    pickle.dump(v, open(fn, 'wb'))
    return v


def shadow_inv(D, J):
    fn = os.path.join(HERE, 'cache', 'bh_shadow.pkl')
    db = pickle.load(open(fn, 'rb')) if os.path.exists(fn) else {}
    key = (round(float(D), 6), J)
    if key not in db:
        import l6_spin_action as L6
        b, logK, resid = L6.shadow(float(D), J)
        db[key] = dict(b=b, logK=logK, resid=resid)
        pickle.dump(db, open(fn, 'wb'))
    return np.linalg.inv(db[key]['b']), db[key]


def helmats(rec, D, J, Kt=None):
    """Mt[nu] (N x 2 x 2) in the helicity basis at Delta = D, cumulative to order Kt."""
    binv, _ = shadow_inv(D, J)
    Ti = HB.basis(J, D)['Ti']
    I = np.array([[EQ.value(rec, 'EE4', a, c, Kt)[0][:N] for c in range(2)] for a in range(2)])      # (a, c, nu)
    M = np.einsum('acn,bc->nab', I, binv)
    return np.einsum('ia,nab,jb->nij', Ti, M, Ti)


def gamma(st, rf, D, J, Kt):
    nus = st['nus'][:N]
    m, eta = BQ.mass(D, J); mref, etaref = BQ.mass(EQ.DREF, J)
    R = BQ.Rj(nus)
    Tst = tcache(nus, eta, R, 'D%g_J%d' % (D, J)); Tref = tcache(nus, etaref, R, 'D%g_J%d' % (EQ.DREF, J))
    NH = sum(R.values())
    cov = sum(R[j]*(nus - 1j*j)**2 for j in R)
    Ms, Mr = helmats(st, D, J, Kt), helmats(rf, EQ.DREF, J)
    out = dict(nus=nus, m=m, NH=NH, cov=cov)
    cs = {}
    for i, h in enumerate(('++', '+-')):
        tot = Tref*Ms[:, i, i]/Mr[:, i, i]
        cc = Tst[0]/tot[0]; cs[h] = cc
        Gex = m*m*(cc*tot - Tst) + REFCORR*GAM[h]*(m/mref)**2*cov
        GK = Gex - GAM[h]*cov
        GK = GK - (GK[0]/NH[0])*NH
        out[h] = dict(GK=GK, Gex=Gex, ratio=Ms[:, i, i]/Mr[:, i, i])
    off = m*m*np.sqrt(cs['++']*cs['+-'])*Tref*0.5*(Ms[:, 0, 1] + Ms[:, 1, 0])/np.sqrt(Mr[:, 0, 0]*Mr[:, 1, 1])
    out['off'] = dict(GK=off, kappa=0.5*(Ms[:, 0, 1] + Ms[:, 1, 0])/np.sqrt(Ms[:, 0, 0]*Ms[:, 1, 1]),
                      asym=0.5*(Ms[:, 0, 1] - Ms[:, 1, 0])/np.sqrt(Ms[:, 0, 0]*Ms[:, 1, 1]))
    return out


if __name__ == '__main__':
    D, J, K, Dg = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    hs = sys.argv[5].split(',') if len(sys.argv) > 5 else ['+-', '++']
    st = pickle.load(open(B4.fn(D, J, K), 'rb')); rf = pickle.load(open(B4.fn(EQ.DREF, J, 2), 'rb'))
    W = BQ.weights()
    lines = []
    def pr(s=''):
        print(s, flush=True); lines.append(s)
    rs = {Kt: gamma(st, rf, D, J, Kt) for Kt in range(0, K + 1)}
    r = rs[K]; m = r['m']
    Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0)); u = (Lam/m)**2
    pr(f'# bh_qan4: EE k=4, Delta={D:.4f} J={J} m={m:.3f} = {m/Lam:.4f} Lam(Dg={Dg}), K={K}, REFCORR={REFCORR}')
    B = HB.basis(J, D)
    pr(f'  helicity basis (order-0 Laplace, x >= 0.15): sv3/sv1 = {B["sv3"]:.1e}, max order-0 offdiag {B["offdiag0"]:.1e}, flatness {B["flat"][0]:.2e} / {B["flat"][1]:.2e}')
    NUP = (0.02, 0.5, 1.02, 1.5, 2.02, 2.5, 3.02, 4.02, 5.02, 6.02, 8.02, 10.02, 13.02, 16.02, 20.02, 30.02, 40.02)
    idx = [int(np.argmin(np.abs(r['nus'] - v))) for v in NUP]
    res = {}
    for h in hs:
        pr(f'  [{h}]  nu: Gamma_K/sum|R_j| (K) [a^K term]  |  Gamma_exact/sum|R_j|  |  HYB gamma_h(nu)^2-scale check: Re Gex/Re cov')
        AR = sum(np.abs(v) for v in BQ.Rj(r['nus']).values())
        for i in idx:
            g = r[h]['GK'][i]; gp = rs[K - 1][h]['GK'][i]
            pr(f'     {r["nus"][i]:5.2f}: {g.real/AR[i]:+10.3f} [{(g - gp).real/AR[i]:+8.3f}]  | {r[h]["Gex"][i].real/AR[i]:+10.3f}'
               f'  | {r[h]["Gex"][i].real/r["cov"][i].real if abs(r["cov"][i].real) > 0 else 0:+8.3f}')
        P = {Kt: GQ.project(dict(nus=r['nus'], GK=rs[Kt][h]['GK']), W) for Kt in range(K + 1)}
        P18 = GQ.project(dict(nus=r['nus'], GK=r[h]['GK']), W, 18.0)
        res[h] = dict(P=P)
        pr(f'  [{h}] projections (nu < 41; [nu < 18]): s_d +- trunc | s_v +- trunc | s_d by order k = 0..K | FF/N (HYB): K=0 -> measured')
        for g in W:
            sd, sv = P[K][g]; ed = abs(sd - P[K - 1][g][0]); ev = abs(sv - P[K - 1][g][1])
            f = GQ.ffn(h, J, u, Dg, sd, sv)
            pr(f'     {g:3s}: s_d = {sd:+9.4f} +- {ed:.4f} [{P18[g][0]:+8.4f}] | s_v = {sv:+9.3f} +- {ev:.3f} | '
               + ' '.join(f'{P[k][g][0] - (P[k-1][g][0] if k else 0):+.4f}' for k in range(K + 1)) + ' | '
               + '  '.join(f'd={dl}: {f[dl][0]:+.4f} -> {f[dl][1]:+.4f} (+- {f[dl][2]*ed:.4f})' for dl in f))
            res[h][g] = dict(sd=sd, ed=ed, sv=sv, ev=ev, ff=f)
        # covariant decomposition (bg_covfit): Gamma_K = dgamma cov + kappa N^_4 + rem, fit nu < 12, weights 1/sum|R_j|
        w = 1/AR; sel = r['nus'] < 12.0
        A = np.array([r['cov'].real[sel]*w[sel], r['NH'].real[sel]*w[sel]]).T
        co, *_ = np.linalg.lstsq(A, r[h]['GK'].real[sel]*w[sel], rcond=None)
        rem = r[h]['GK'] - co[0]*r['cov'] - co[1]*r['NH']
        Pr = GQ.project(dict(nus=r['nus'], GK=rem), W)
        pr(f'  [{h}] covariant part dgamma = {co[0]:+.4f} (= {co[0]/GAM[h]*100:+.2f} % of gamma_h); remainder s_d: '
           + ' '.join(f'{g}:{Pr[g][0]:+.4f}' for g in W) + f'; max |rem|/sum|R_j| nu<2.5 / 2.5-18: '
           f'{np.max(np.abs(rem*w)[r["nus"] < 2.5]):.2f} / {np.max(np.abs(rem*w)[(r["nus"] >= 2.5) & (r["nus"] < 18)]):.2f}')
        res[h]['cov'] = dict(dgamma=co[0], rem={g: Pr[g][0] for g in W})
    # helicity mixing
    k = r['off']['kappa']; AR = sum(np.abs(v) for v in BQ.Rj(r['nus']).values())
    pr('  mixing: m^2 kappa(nu) (state, symmetric / antisymmetric part), Gamma_off/sum|R_j| at nu = ' + ' '.join(f'{r["nus"][i]:.2f}' for i in idx))
    pr('     m^2 kappa  : ' + ' '.join(f'{(m*m*k[i]).real:+.2f}' for i in idx))
    pr('     m^2 asym   : ' + ' '.join(f'{(m*m*r["off"]["asym"][i]).real:+.2f}' for i in idx))
    pr('     Gamma_off/R: ' + ' '.join(f'{(r["off"]["GK"][i]/AR[i]).real:+.2f}' for i in idx))
    Po = {Kt: GQ.project(dict(nus=r['nus'], GK=rs[Kt]['off']['GK']), W) for Kt in (K - 1, K)}
    pr('     s_off (derivative projection of Re Gamma_off) per gauge: ' + ' '.join(f'{g}:{Po[K][g][0]:+.4f}(+-{abs(Po[K][g][0] - Po[K-1][g][0]):.4f})' for g in W))
    res['off'] = {g: Po[K][g] for g in W}; res['off_trunc'] = {g: Po[K - 1][g] for g in W}
    open(os.path.join(HERE, 'cache', 'bh_qan4.txt'), 'a').write('\n'.join(lines) + '\n')
    pickle.dump(dict(rs=rs, res=res, D=D, J=J, Dg=Dg, m=m, u=u), open(os.path.join(HERE, 'cache', 'bh_qan4_D%g_J%d_Dg%d.pkl' % (D, J, Dg)), 'wb'))
