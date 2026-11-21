"""
bg_qan (BELOW_FOLD2.md §4): Result Q chain of bf_qan (read-only import: mass, Rj, T, weights) on the analytic-radial EE k = 4 images of
bg_eeq.  Differences from bf_qan: the image is exact order by order (no QMC); the normalisation N(nu) = I_ref/T_ref needs no fit
(total_st/(2c) = T_ref(nu) * I_st/I_ref); the projections run over the whole de-shift grid (nu < 41).
  Gamma_exact = m^2 (c T_ref I_st/I_ref - T_st),  c: Gamma_exact(nu_min) = 0;   Gamma_K = Gamma_exact - Gamma_HYB - const N^_4 (Gamma_K(nu_min)=0),
  Gamma_HYB = gamma sum_j R_j (nu - ij)^2, gamma = (c4 + 1/2)/2 = 0.75 ('++'), 20.75 ('+-').
  s_d = 0.04 sum chi_d Re Gamma_K;  s_v likewise (a_h-convention dependent);  truncation error = |s(K) - s(K-1)| (size of the a^K term).
  FF/N at the measured state (Stage-6 HYB model, stage6_model.hyb):  FF -> FF - a 2u^3 s_d - c (2u^3/Lam^2) s_v.
  effective K: s_d / s_a, s_d / s_b with Stage 7's derivative stress coefficients (structures a, b).
Usage: python3 bg_qan.py Delta J K Dg h   -> appends to cache/bg_qan.txt, pickle cache/bg_qan_D{Delta}_J{J}.pkl
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

GAM = {'++': 0.75, '+-': 20.75}


def tcache(nus, eta, R, key):
    fn = os.path.join(HERE, 'cache', 'bg_qan_T_%s.pkl' % key)
    if os.path.exists(fn):
        d = pickle.load(open(fn, 'rb'))
        if len(d) == len(nus):
            return d
    v = BQ.T(nus, eta, R)
    pickle.dump(v, open(fn, 'wb'))
    return v


def gamma_K(db, h='++', kind='EE4', a=0, c=0, Kt=None, ref_kind=None):
    st, rf = db['st'], db['ref']
    D, J = st['Delta'], st['J']
    n = 1024
    nus = st['nus'][:n]
    m, eta = BQ.mass(D, J); mref, etaref = BQ.mass(EQ.DREF, J)
    R = BQ.Rj(nus)
    Tst = tcache(nus, eta, R, 'D%g_J%d' % (D, J)); Tref = tcache(nus, etaref, R, 'D%g_J%d' % (EQ.DREF, J))
    NH = sum(R.values())
    HYB = GAM[h]*sum(R[j]*(nus - 1j*j)**2 for j in R)
    Ist, pst = EQ.value(st, kind, a, c, Kt)
    Iref, _ = EQ.value(rf, ref_kind or kind, a, c)
    Ist, Iref = Ist[:n], Iref[:n]
    tot = Tref*Ist/Iref
    cc = Tst[0]/tot[0]
    Gex = m*m*(cc*tot - Tst)
    GK = Gex - HYB
    GK = GK - (GK[0]/NH[0])*NH
    return dict(nus=nus, GK=GK, Gex=Gex, m=m, NH=NH, Tst=Tst, Tref=Tref, ratio=Ist/Iref)


def project(r, W, numax=41.0):
    out = {}
    for g, (nu, cv, cd) in W.items():
        k = min(len(r['nus']), len(cd))
        sel = r['nus'][:k] < numax
        sd = 0.04*np.sum(cd[:k][sel]*r['GK'][:k][sel].real)
        sv = 0.04*np.sum(cv[:k][sel]*r['GK'][:k][sel].real)
        out[g] = (sd, sv)
    return out


STRESS = {'H0': (-0.90, +0.37), 'H2': (-0.83, +0.04), 'H3': (-0.90, +0.38), 'H4': (-0.89, -0.01), 'SB': (0.0, -1.17), 'XP': (-0.58, -0.17)}


def ffn(h, J, u, Dg, sd, sv):
    """FF/N (Stage-6 HYB) at the state, K = 0 and with the measured (s_d, s_v), for delta = 0.05, 0.2."""
    import stage6_model as M, common_q as CQ
    from stage6_scan import BMAX
    out = {}
    for delta in (0.05, 0.2):
        name = M.FUNCS[delta]; F, c2, c3 = M.load(name); a, c = F['a'], F['c']
        Lam = M.lam(Dg); tab = M.Tables(name, Lam, BMAX[Dg])
        Js, B2, B3 = CQ.parts(h, u, max(J, 4), [c2], [c3]); k = list(Js).index(J)
        FF, N, _, _ = M.hyb(F, h, u, float(J), B2[k, 0], B3[k, 0], tab, Lam)
        dFF = -a*2*u**3*sd - c*(2*u**3/Lam**2)*sv
        out[delta] = (float(FF/N), float((FF + dFF)/N), float(a*2*u**3/N), float(c*2*u**3/Lam**2/N))
    return out


if __name__ == '__main__':
    D, J, K = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    Dg = int(sys.argv[4]); h = sys.argv[5] if len(sys.argv) > 5 else '++'
    db = pickle.load(open(EQ.fn(D, J, K), 'rb'))
    W = BQ.weights()
    lines = []
    def pr(s=''):
        print(s, flush=True); lines.append(s)
    rs = {Kt: gamma_K(db, h, Kt=Kt) for Kt in range(0, K + 1)}
    r = rs[K]
    u = (np.sqrt((Dg - 2.0)*(Dg - 1.0))/r['m'])**2
    pr(f'# bg_qan: EE k=4 analytic-radial image, Delta={D:.4f} J={J} ({h}) m={r["m"]:.3f} = {r["m"]/np.sqrt((Dg-2.)*(Dg-1.)):.4f} Lam(Dg={Dg}), K={K}')
    pr('  nu     Gamma_K(K)   [a^K-term]   Gamma_K(K-1)  |  Gamma_exact')
    for i, v in enumerate(r['nus']):
        if v in (0.02, 0.5, 1.02, 1.5, 2.02, 2.5, 3.02, 4.02, 5.02, 6.02, 8.02, 10.02, 13.02, 16.02, 20.02, 25.02, 30.02, 40.02) or abs(v - round(v)) < 1e-9:
            if v in (0.02, 0.5, 1.02, 1.5, 2.02, 2.5, 3.02, 4.02, 5.02, 6.02, 8.02, 10.02, 13.02, 16.02, 20.02, 25.02, 30.02, 40.02):
                pr(f'  {v:5.2f}  {r["GK"][i].real:+10.3f}  [{(r["GK"][i] - rs[K-1]["GK"][i]).real:+9.3f}]  {rs[K-1]["GK"][i].real:+10.3f}'
                   f'  |  {r["Gex"][i].real:+10.3f}')
    P = {Kt: project(rs[Kt], W) for Kt in range(0, K + 1)}
    P18 = project(r, W, 18.0)
    pr('  projections (whole grid nu < 41; nu < 18 in brackets): s_d +- trunc | s_v +- trunc | K_eff (a), (b) | s_d by order k = 0..K')
    for g in W:
        sd, sv = P[K][g]; ed = abs(sd - P[K - 1][g][0]); ev = abs(sv - P[K - 1][g][1])
        sa, sb = STRESS[g]
        Ka = sd/sa if sa else np.nan
        pr(f'     {g:3s}: s_d = {sd:+9.4f} +- {ed:.4f} [{P18[g][0]:+9.4f}] | s_v = {sv:+10.3f} +- {ev:.3f} | K_eff = {Ka:+8.3f}, {sd/sb:+8.3f}'
           f' | ' + ' '.join(f'{P[k][g][0] - (P[k-1][g][0] if k else 0):+.4f}' for k in range(K + 1)))
    pr('  FF/N (HYB, Stage 6) at this state: delta: K=0 -> with measured Gamma_K, per gauge')
    for g in W:
        sd, sv = P[K][g]
        f = ffn(h, J, u, Dg, sd, sv)
        ed = abs(sd - P[K - 1][g][0])
        pr(f'     {g:3s}: ' + '  '.join(f'd={dl}: {f[dl][0]:+.4f} -> {f[dl][1]:+.4f} (+- {f[dl][2]*ed:.4f})' for dl in f))
    open(os.path.join(HERE, 'cache', 'bg_qan.txt'), 'a').write('\n'.join(lines) + '\n')
    pickle.dump(dict(rs=rs, P=P), open(os.path.join(HERE, 'cache', 'bg_qan_D%g_J%d.pkl' % (D, J)), 'wb'))
