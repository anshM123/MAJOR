"""
po_qan5 (PARITY_ODD_CLOSE.md sect. 2): Result Q chain (bh_qan4 / bg_qan logic, copied) for the PARITY-ODD +- state at odd J: the 1 x 1
image of po_ee_odd (symmetric parity-odd <TTO> structure on both sides; the shadow coefficient is a nu-independent constant and drops out
with the per-state constant c).
  Gamma_exact = m^2 (c T_ref I_st/I_ref - T_st) + gamma (m/m_ref)^2 sum_j R_j (nu - ij)^2   (REFCORR = 1 as bh_qan4),
  Gamma_K = Gamma_exact - gamma sum_j R_j (nu - ij)^2 - const N^_4,  gamma = gamma_{+-} = 20.75 (HYB, helicity-only),
  s_d, s_v per gauge (bf_qan.weights), truncation = |s(K) - s(K-1)|; FF/N at the state (bg_qan.ffn, Stage-6 HYB, h = '+-').
Usage: python3 po_qan5.py Delta J K Dg   (needs cache/po_eeq_D{Delta}_J{J}_K{K}.pkl and po_eeq_D20000_J{J}_K2.pkl)
       -> appends cache/po_qan5.txt, pickle cache/po_qan5_D{Delta}_J{J}_Dg{Dg}.pkl (bh_ffscan-compatible 'res' layout)
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
RQ = os.path.join(HERE, '..', 'resultq'); sys.path.insert(0, RQ)
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import bf_qan as BQ
import bg_eeq as EQ
import bg_qan as GQ
import po_ee_odd as PO

GAM = 20.75
REFCORR = float(os.environ.get('PO_REFCORR', '1'))
N = 1024


def tcache(nus, eta, R, key):
    for fn in (os.path.join(HERE, 'cache', 'bh_qan_T_%s.pkl' % key), os.path.join(HERE, 'cache', 'po_qan_T_%s.pkl' % key)):
        if os.path.exists(fn):
            d = pickle.load(open(fn, 'rb'))
            if len(d) == len(nus):
                return d
    v = BQ.T(nus, eta, R)
    pickle.dump(v, open(os.path.join(HERE, 'cache', 'po_qan_T_%s.pkl' % key), 'wb'))
    return v


def gamma(st, rf, D, J, Kt):
    nus = st['nus'][:N]
    m, eta = BQ.mass(D, J); mref, etaref = BQ.mass(EQ.DREF, J)
    R = BQ.Rj(nus)
    Tst = tcache(nus, eta, R, 'D%g_J%d' % (D, J)); Tref = tcache(nus, etaref, R, 'D%g_J%d' % (EQ.DREF, J))
    NH = sum(R.values())
    cov = sum(R[j]*(nus - 1j*j)**2 for j in R)
    Is = EQ.value(st, 'EE4', 0, 0, Kt)[0][:N]; Ir = EQ.value(rf, 'EE4', 0, 0)[0][:N]
    tot = Tref*Is/Ir
    cc = Tst[0]/tot[0]
    Gex = m*m*(cc*tot - Tst) + REFCORR*GAM*(m/mref)**2*cov
    GK = Gex - GAM*cov
    GK = GK - (GK[0]/NH[0])*NH
    return dict(nus=nus, m=m, NH=NH, cov=cov, GK=GK, Gex=Gex, ratio=Is/Ir)


if __name__ == '__main__':
    D, J, K, Dg = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    st = pickle.load(open(PO.fn(D, J, K), 'rb')); rf = pickle.load(open(PO.fn(EQ.DREF, J, 2), 'rb'))
    W = BQ.weights()
    lines = []

    def pr(s=''):
        print(s, flush=True); lines.append(s)
    rs = {Kt: gamma(st, rf, D, J, Kt) for Kt in range(0, K + 1)}
    r = rs[K]; m = r['m']
    Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0)); u = (Lam/m)**2
    pr(f'# po_qan5: parity-odd +- EE k=4, Delta={D:.4f} J={J} m={m:.3f} = {m/Lam:.4f} Lam(Dg={Dg}), K={K}, REFCORR={REFCORR}')
    NUP = (0.02, 0.5, 1.02, 1.5, 2.02, 2.5, 3.02, 4.02, 5.02, 6.02, 8.02, 10.02, 13.02, 16.02, 20.02, 30.02, 40.02)
    idx = [int(np.argmin(np.abs(r['nus'] - v))) for v in NUP]
    AR = sum(np.abs(v) for v in BQ.Rj(r['nus']).values())
    pr('  nu: Gamma_K/sum|R_j| (K) [a^K term]  |  Gamma_exact/sum|R_j|  |  Re Gex/Re cov')
    for i in idx:
        g = r['GK'][i]; gp = rs[K - 1]['GK'][i]
        pr(f'     {r["nus"][i]:5.2f}: {g.real/AR[i]:+10.3f} [{(g - gp).real/AR[i]:+8.3f}]  | {r["Gex"][i].real/AR[i]:+10.3f}'
           f'  | {r["Gex"][i].real/r["cov"][i].real if abs(r["cov"][i].real) > 0 else 0:+8.3f}')
    P = {Kt: GQ.project(dict(nus=r['nus'], GK=rs[Kt]['GK']), W) for Kt in range(K + 1)}
    P18 = GQ.project(dict(nus=r['nus'], GK=r['GK']), W, 18.0)
    res = {'+-': {}}
    pr('  projections (nu < 41; [nu < 18]): s_d +- trunc | s_v +- trunc | s_d by order k = 0..K | FF/N (HYB): K=0 -> measured')
    for g in W:
        sd, sv = P[K][g]; ed = abs(sd - P[K - 1][g][0]); ev = abs(sv - P[K - 1][g][1])
        f = GQ.ffn('+-', J, u, Dg, sd, sv)
        pr(f'     {g:3s}: s_d = {sd:+9.4f} +- {ed:.4f} [{P18[g][0]:+8.4f}] | s_v = {sv:+9.3f} +- {ev:.3f} | '
           + ' '.join(f'{P[k][g][0] - (P[k-1][g][0] if k else 0):+.4f}' for k in range(K + 1)) + ' | '
           + '  '.join(f'd={dl}: {f[dl][0]:+.4f} -> {f[dl][1]:+.4f} (+- {f[dl][2]*ed:.4f})' for dl in f))
        res['+-'][g] = dict(sd=sd, ed=ed, sv=sv, ev=ev, ff=f)
    w = 1/AR; sel = r['nus'] < 12.0
    A = np.array([r['cov'].real[sel]*w[sel], r['NH'].real[sel]*w[sel]]).T
    co, *_ = np.linalg.lstsq(A, r['GK'].real[sel]*w[sel], rcond=None)
    rem = r['GK'] - co[0]*r['cov'] - co[1]*r['NH']
    Pr = GQ.project(dict(nus=r['nus'], GK=rem), W)
    pr(f'  covariant part dgamma = {co[0]:+.4f} (= {co[0]/GAM*100:+.2f} % of gamma_+-); remainder s_d: '
       + ' '.join(f'{g}:{Pr[g][0]:+.4f}' for g in W))
    # cross-check against the universal helicity-odd correction of BELOW_FOLD3 (++ J = 0/2/4 and +- J = 4 images, sigma_+- = -1)
    try:
        ref4 = pickle.load(open(os.path.join(HERE, 'cache', 'bh_qan4_D%g_J4_Dg%d.pkl' % ({300: 304.0297506031432, 1000: 1004.0089775159122}[Dg], Dg)), 'rb'))
        g4 = ref4['rs'][max(ref4['rs'])]['+-']['GK']
        AR4 = AR
        pr('  cross-check vs +- J = 4 (BELOW_FOLD3, same Dg), Gamma_K/sum|R_j| at nu = ' + ' '.join(f'{r["nus"][i]:.2f}' for i in idx[:12]))
        pr('     J=5 odd : ' + ' '.join(f'{(r["GK"][i]/AR[i]).real:+.2f}' for i in idx[:12]))
        pr('     J=4 even: ' + ' '.join(f'{(g4[i]/AR4[i]).real:+.2f}' for i in idx[:12]))
        pr('     s_d (J=4 +-): ' + ' '.join(f'{gg}:{ref4["res"]["+-"][gg]["sd"]:+.4f}' for gg in W))
    except Exception as e:
        pr(f'  (no J = 4 cross-check: {e})')
    open(os.path.join(HERE, 'cache', 'po_qan5.txt'), 'a').write('\n'.join(lines) + '\n')
    pickle.dump(dict(rs=rs, res=res, D=D, J=J, Dg=Dg, m=m, u=u), open(os.path.join(HERE, 'cache', 'po_qan5_D%g_J%d_Dg%d.pkl' % (D, J, Dg)), 'wb'))
