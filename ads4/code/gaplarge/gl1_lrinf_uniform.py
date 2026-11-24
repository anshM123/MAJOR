"""
gaplarge/gl1_lrinf_uniform.py (GAP_LARGESPIN.md §1): Lemma LR-inf for EVERY Delta_g >= 1000, certified.

Scaling lemma ([P], GAP_LARGESPIN §1.2): with the envelope R^ (gl_common), Lam^4 x N(beta; Lam) is non-increasing in Lam, where
N = (sqrt2/pi) sum_k w_k [Ta + Tb + Tc1^ + Tc2] >= 2|E2| + |E3|.  Hence for every Lam >= Lam0 = Lam(1000) and every beta:
     r(beta; Lam)      <= (Lam0/Lam)^2 * r^0(beta),        r^0  := Lam0^2 N(beta; Lam0)/(A1min M1(beta))
     rtilde(beta; Lam) <= (Lam0/Lam)   * rt^0(beta),       rt^0 := Lam0^2 N(beta; Lam0)/(A1min M1 + A2' Lam0^2/(2 bhat0^3))
     Lam^2 F_LR(beta, u, sigma; Lam) >= min( A1min M1(beta), Lam0^2 F^0_cert(beta) )   (u <= 1, sigma = +-1)
All three right-hand sides are evaluated here cell-wise in Arb at Lam0 (cells of lrt3: 830 on [0.1, 24], 80 on (24, 64]; tail
beta >= 64).  Optional argv: extra gaps at which N is ALSO evaluated directly (check of the scaling bound, e.g. 2000 5000).
Output: cache/gl1_lrinf_uniform.{json,txt}
"""
import os, sys, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from flint import arb
import gl_common as G
from gl_common import A, up, lo, W, LAM0, A1MIN, A2P

lines = []
def pr(s=''):
    print(s, flush=True); lines.append(s)


def run_cells(LAM, cl, TAB):
    rows = []
    for (a, b, m1) in cl:
        assert m1 > 0
        e2, p2 = G.E_up(2, a, b, TAB, LAM); e3, _ = G.E_up(3, a, b, TAB, LAM)
        num = 2*e2 + e3
        main_ = A(A1MIN)*m1/LAM**2
        sb = 2*(A(b)/2).sinh()                     # bhat/Lam at the cell's right end (flat term decreasing)
        flat = A(A2P)/(2*(LAM*sb)**3)
        r = up(num/main_); rt = up(num/(main_ + flat)); Fc = lo(main_ + flat - num)
        rows.append(dict(a=a, b=b, m1=m1, num=num, r=r, rt=rt, Fc=Fc, L2Fc=lo(LAM**2*(main_ + flat - num)),
                         A1M1=lo(A(A1MIN)*m1), parts2={k: float(up(v).mid()) for k, v in p2.items()}))
    return rows


def sup_from(rows, key, bp, tail):
    v = tail
    for R_ in rows:
        if R_['a'] >= bp:
            v = v if v.upper() >= R_[key].upper() else R_[key]
    return v


def beta0(rows, key, thr):
    """smallest cell start beta0 with every cell from beta0 on having key < thr (Arb-decided; tail checked by caller)."""
    b0 = None
    for R_ in reversed(rows):
        if R_[key] < thr:
            b0 = R_['a']
        else:
            break
    return b0


def main():
    pr(__doc__.split('Output')[0])
    t0 = time.time()
    for k in (2, 3):
        pr(f'envelope k = {k}: s1 = {W[k].s1.str(8)}, V = sup s^4 R(s) <= {W[k].V.str(8)}  '
           f'(c3 = |psi\'\'\'(0)| = {float(abs(W[k].d0[3])):.6g}; R^ = V/s^4 below s1, = R_inf above)')
    cl = G.m1_cells()
    TAB0 = {k: G.build_I(W[k], LAM0) for k in (2, 3)}
    pr(f'I-tables at Lam0 built ({time.time()-t0:.0f}s)')
    rows = run_cells(LAM0, cl, TAB0)
    tail0, ell = G.tail_majorant(TAB0, LAM0)
    pr(f'Lam0 = {LAM0.str(12)}; tail: r^0(beta) <= {tail0.str(5)} for all beta >= 64  ({time.time()-t0:.0f}s)')

    # comparison with LRINF_TIGHT's R-based cell values (cost of the envelope)
    old = json.load(open(os.path.join(G.LRT_CACHE, 'lrt3_rcert.json')))['rows']
    oldr = {(Fr(o['ba']), Fr(o['bb'])): o for o in old}
    pr('\ncell values at Lam0 with the envelope R^ (this note) vs R (LRINF_TIGHT): r | rtilde | F_cert')
    for R_ in rows:
        a = R_['a']
        if a in (Fr('0.1'), Fr('0.15'), Fr('0.2'), Fr('0.245'), Fr('0.3'), Fr('0.35'), Fr('0.365'), Fr('0.37'), Fr('0.4'), Fr('0.5'),
                 Fr('0.7'), Fr('1'), Fr('1.5'), Fr('2'), Fr('3'), Fr('5'), Fr('10'), Fr('20'), Fr('40'), Fr('63.5')):
            o = oldr[(R_['a'], R_['b'])]
            pr(f"  [{float(a):6.3f},{float(R_['b']):6.3f}]  r^0 <= {float(R_['r'].mid()):.4e} (R: {o['r']:.4e})   rt^0 <= "
               f"{float(R_['rt'].mid()):.4e} (R: {o['rt']:.4e})   F^0 >= {float(R_['Fc'].mid()):.4e} (R: {o['Fcert']:.4e})   "
               f"[k=2 Tc1 {R_['parts2']['Tc1']:.2e} Tc2 {R_['parts2']['Tc2']:.2e}]")

    # uniform thresholds
    pr('\nsup over beta >= beta\' of r^0 and rt^0 (cells + tail):')
    for bp in ['0.1', '0.15', '0.2', '0.245', '0.25', '0.3', '0.35', '0.37', '0.4', '0.44', '0.5', '0.6', '0.8', '1', '1.5', '2', '3', '5']:
        bp = Fr(bp)
        pr(f"  beta >= {float(bp):5.3f}: sup r^0 <= {float(sup_from(rows, 'r', bp, tail0).mid()):.4e}   "
           f"sup rt^0 <= {float(sup_from(rows, 'rt', bp, tail0).mid()):.4e}")
    out = dict(Lam0=float(LAM0.mid()), tail=float(tail0.mid()), V={k: float(W[k].V.mid()) for k in (2, 3)},
               s1={k: float(W[k].s1.mid()) for k in (2, 3)})
    assert tail0 < 1
    for key in ('r', 'rt'):
        for thr in (1, Fr(1, 2)):
            b0 = beta0(rows, key, A(thr))
            out[f'beta0_{key}_{float(thr)}'] = str(b0)
            pr(f'uniform (all Dg >= 1000): smallest certified beta0 with {key}^0 < {float(thr)} on [beta0, oo): {float(b0)}')

    # Lam-improved thresholds: r(beta; Lam) <= (Lam0/Lam)^2 r^0 < 1  <=  r^0 < (Lam/Lam0)^2 ; rt <= (Lam0/Lam) rt^0
    pr('\nLam-improved certified thresholds beta0(Dg) (r < 1 form: r^0 < (Lam/Lam0)^2;  F_LR > 0 form: rt^0 < Lam/Lam0 or '
       'r^0 < (Lam/Lam0)^2), beta-cells stop at 0.1:')
    impr = {}
    for Dg in (1000, 1500, 2000, 3000, 5000, 10000, 20000, 50000):
        L = G.lam_of(Dg); q = L/LAM0
        b_r = beta0(rows, 'r', q**2)
        b_t = None
        for R_ in reversed(rows):
            if R_['rt'] < q or R_['r'] < q**2:
                b_t = R_['a']
            else:
                break
        impr[Dg] = dict(Lam=float(L.mid()), beta0_r=str(b_r), beta0_rt=str(b_t),
                        s0_r=float((2*(A(b_r)/2).sinh()).mid()), s0_rt=float((2*(A(b_t)/2).sinh()).mid()))
        pr(f'  Dg = {Dg:6d} (Lam/Lam0 = {float(q.mid()):.4f}): beta0(r<1) = {float(b_r):.4f} (s0 = bhat0/Lam = {impr[Dg]["s0_r"]:.5f}, '
           f'bhat0 = {impr[Dg]["s0_r"]*float(L.mid()):.1f});  beta0(F_LR>0) = {float(b_t):.4f} (s0 = {impr[Dg]["s0_rt"]:.5f})')
    out['improved'] = impr

    # uniform gap-unit floor
    pr('\nuniform gap-unit floor  Lam^2 F_LR >= f_unif(beta) := min(A1min M1, Lam0^2 F^0_cert)  (all Dg >= 1000, u <= 1, sigma):')
    for bp in ['0.245', '0.3', '0.37', '0.4', '0.5', '0.7', '1', '1.5', '2', '3', '5', '7.5', '10', '20']:
        bp = Fr(bp)
        R_ = next(R for R in rows if R['a'] == bp)
        fu = R_['L2Fc'] if R_['L2Fc'] < R_['A1M1'] else R_['A1M1']
        pr(f"  cell [{float(R_['a']):.3f}, {float(R_['b']):.3f}]: A1min M1 >= {float(R_['A1M1'].mid()):.5f};  Lam0^2 F^0 >= "
           f"{float(R_['L2Fc'].mid()):.5f};  f_unif >= {float(fu.mid()):.5f}   (r^0 <= {float(R_['r'].mid()):.4f})")

    rows_out = [dict(ba=str(R_['a']), bb=str(R_['b']), M1lo=float(R_['m1'].mid()), r0=float(R_['r'].mid()),
                     rt0=float(R_['rt'].mid()), F0=float(R_['Fc'].mid()), L2F0=float(R_['L2Fc'].mid()),
                     A1M1=float(R_['A1M1'].mid()), L2N0=float(up(LAM0**2*R_['num']).mid())) for R_ in rows]
    out['rows'] = rows_out

    # optional direct evaluation at other gaps: check N(beta; Lam) <= (Lam0/Lam)^4 N(beta; Lam0) cell-wise
    for Dg in [int(x) for x in sys.argv[1:]]:
        L = G.lam_of(Dg)
        TAB = {k: G.build_I(W[k], L) for k in (2, 3)}
        rr = run_cells(L, cl, TAB)
        tl, _ = G.tail_majorant(TAB, L)
        worst = 0.0; viol = 0
        for R0, R1 in zip(rows, rr):
            bound = (LAM0/L)**2*R0['r']
            if not (R1['r'] <= bound):
                viol += 1
            worst = max(worst, float((R1['r']/bound).mid()))
        b_dir = beta0(rr, 'r', arb(1)); bt_dir = beta0(rr, 'rt', arb(1))
        pr(f'\ndirect evaluation at Dg = {Dg} (Lam = {L.str(10)}): r(Lam) / [(Lam0/Lam)^2 r^0] max over cells = {worst:.4f} '
           f'(violations: {viol}); tail {tl.str(4)} vs scaled {float(((LAM0/L)**2*tail0).mid()):.3e};  direct beta0(r<1) = '
           f'{float(b_dir)}, beta0(rt<1) = {float(bt_dir)}')
        for bp in ['0.1', '0.2', '0.3', '0.37', '0.5', '1', '2']:
            bp = Fr(bp)
            pr(f"   beta >= {float(bp):4.2f}: sup r direct {float(sup_from(rr, 'r', bp, tl).mid()):.4e}  vs (Lam0/Lam)^2 sup r^0 "
               f"{float(((LAM0/L)**2*sup_from(rows, 'r', bp, tail0)).mid()):.4e};  sup rt direct {float(sup_from(rr, 'rt', bp, tl).mid()):.4e}")
        out[f'direct_{Dg}'] = dict(worst_ratio=worst, violations=viol, beta0_r=str(b_dir), beta0_rt=str(bt_dir))
    json.dump(out, open(os.path.join(G.CACHE, 'gl1_lrinf_uniform.json'), 'w'), indent=0)
    open(os.path.join(G.CACHE, 'gl1_lrinf_uniform.txt'), 'w').write('\n'.join(lines) + '\n')
    pr(f'done in {time.time()-t0:.0f}s')


if __name__ == '__main__':
    main()
