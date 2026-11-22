"""
rs2_trunc (HRS_CLOSE.md §4-5): with H_RS reduced to the engine's own truncation (HRS_CLOSE §2), the remaining error of the binding
state's s_d is the omitted a^{K+1} term.  This script reads the K = 5 image (rs2_engine5) and reports, per nu window:
  * the per-harmonic a^4 and a^5 terms of Re Gamma_K / sum_j |R_j(nu)| (c2 per harmonic, Result Q's unit);
  * the window tolerances: the uniform per-harmonic |dG| that would use up the delta = 0.2 / 0.05 margins (H0 weight
    0.04 |chi_d| sum|R_j|, margins s_fail - s_d = 1.78 - 0.51 and 0.58 - 0.51, as HRS_CHECK §4);
  * s_d (all six gauges) at K = 4 and K = 5, and FF/N (bg_qan.ffn, Stage-6 HYB) with both.
Read-only imports: bg_qan (gamma_K, project, ffn), bf_qan (Rj, weights), bg_eeq.
Usage: python3 rs2_trunc.py Dg [full]   -> cache/rs2_trunc_D{Dg}[full].txt
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'resultq'))
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import bg_qan as BG, bf_qan as BQ
import rs2_engine5 as E5

SFAIL = {0.2: 1.78, 0.05: 0.58}
WIN = ((0, 1), (1, 2), (2, 3), (3, 4), (4, 6), (6, 8), (8, 10), (10, 18), (18, 41))


def main(Dg, full=False):
    lines = []
    def pr(s=''):
        print(s, flush=True); lines.append(s)
    db = E5.get(float(Dg), 0, full)
    old = pickle.load(open(BG.EQ.fn(float(Dg), 0, 4), 'rb'))
    rs = {Kt: BG.gamma_K(db, '++', Kt=Kt) for Kt in range(0, 6)}
    r4o = BG.gamma_K(old, '++', Kt=4)
    nus = rs[5]['nus']; R = BQ.Rj(nus); SR = sum(np.abs(v) for v in R.values())
    W = BQ.weights()
    chk = np.max(np.abs(rs[4]['GK'] - r4o['GK'])/np.maximum(SR, 1e-300))
    pr(f'# rs2_trunc: Delta_g = {Dg}, ++ J = 0 at the gap (Delta = {Dg}), EE k = 4; K = 5 state image (rs2_eeq), reference from bg_eeq K4')
    pr(f'  check: K = 4 truncation of the new image vs the production bg_eeq K = 4 image: max |dGamma_K|/sum|R_j| = {chk:.2e}')
    a4 = ((rs[4]['GK'] - rs[3]['GK'])/SR).real
    a5 = ((rs[5]['GK'] - rs[4]['GK'])/SR).real
    g5 = (rs[5]['GK']/SR).real
    _, _, cd = W['H0']
    k = min(len(nus), len(cd))
    wgt = 0.04*np.abs(np.asarray(cd[:k]))*SR[:k]
    sd = {Kt: BG.project(rs[Kt], W) for Kt in range(0, 6)}
    s0 = sd[5]['H0'][0]
    pr('  window   | weight share | tolerance (uniform |dG| per harmonic) d=0.2 / d=0.05 | max|Gamma_K| | max|a^4 term| | max|a^5 term| | a^5 term in s_d (H0)')
    tot = wgt.sum()
    for (lo, hi) in WIN:
        m = (nus[:k] >= lo) & (nus[:k] < hi)
        if not m.any():
            continue
        w = wgt[m].sum()
        t02 = (SFAIL[0.2] - s0)/w; t005 = (SFAIL[0.05] - s0)/w
        ds5 = 0.04*np.sum(np.asarray(cd[:k])[m]*SR[:k][m]*a5[:k][m])
        pr(f'  [{lo:2d},{hi:2d})  | {100*w/tot:5.1f} %      | {t02:9.1f} / {t005:7.2f}                           | {np.abs(g5[:k][m]).max():9.2f}  '
           f'| {np.abs(a4[:k][m]).max():9.4f}    | {np.abs(a5[:k][m]).max():9.5f}    | {ds5:+.5f}')
    pr('  s_d per gauge: K = 3 / K = 4 / K = 5 | a^4 term | a^5 term')
    for g in W:
        pr(f'     {g:3s}: {sd[3][g][0]:+.5f} / {sd[4][g][0]:+.5f} / {sd[5][g][0]:+.5f} | {sd[4][g][0] - sd[3][g][0]:+.5f} | {sd[5][g][0] - sd[4][g][0]:+.6f}')
    u = (np.sqrt((Dg - 2.0)*(Dg - 1.0))/rs[5]['m'])**2
    pr('  FF/N (Stage-6 HYB, bg_qan.ffn) with s_d, s_v at K = 4 and K = 5 (H0 / SB)')
    for g in ('H0', 'SB'):
        f4 = BG.ffn('++', 0, u, int(Dg), *sd[4][g]); f5 = BG.ffn('++', 0, u, int(Dg), *sd[5][g])
        pr(f'     {g}: ' + '  '.join(f'd={dl}: K4 {f4[dl][1]:+.5f}  K5 {f5[dl][1]:+.5f}  (diff {f5[dl][1] - f4[dl][1]:+.6f})' for dl in f4))
    txt = '\n'.join(lines)
    open(os.path.join(HERE, 'cache', f'rs2_trunc_D{Dg}{"full" if full else ""}.txt'), 'w').write(txt + '\n')
    pickle.dump(dict(nus=nus, a4=a4, a5=a5, g5=g5, sd=sd), open(os.path.join(HERE, 'cache', f'rs2_trunc_D{Dg}{"full" if full else ""}.pkl'), 'wb'))


if __name__ == '__main__':
    main(int(sys.argv[1]), len(sys.argv) > 2 and sys.argv[2] == 'full')
