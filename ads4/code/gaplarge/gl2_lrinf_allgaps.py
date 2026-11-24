"""
gaplarge/gl2_lrinf_allgaps.py (GAP_LARGESPIN.md §1.3-1.5): Lemma LR-inf CERTIFIED for every Delta_g >= 1000 with LRINF_TIGHT's own
thresholds (beta0 = 0.37 for r < 1, 0.245 for F_LR > 0), by
  Segment A (Lam0 <= Lam <= Lam_n): an adaptive geometric grid Lam_0 < Lam_1 < ... < Lam_n (exact rationals; Lam_0 a rational lower
     bound of Lam(1000), Lam_n >= Lam(2000)).  Every piece of LRINF_TIGHT's bound N(beta; Lam) (with R itself) is NON-INCREASING in Lam
     at fixed beta [P] (Ta, Tb, Tc2: positive coefficients x negative powers of Lam; Tc1: R(Lam t) decreasing in Lam since both
     branches of R are decreasing in x).  Hence for Lam in [Lam_i, Lam_{i+1}]:
        r(beta; Lam)  <= Lam_{i+1}^2 N(beta; Lam_i)/(A1min M1)
        rt(beta; Lam) <= Lam_{i+1}^2 N(beta; Lam_i)/(A1min M1 + A2'/(2 Lam_{i+1} s^3)),     s = bhat/Lam = 2 sinh(beta/2)
        Lam^2 F_LR    >= A1min M1 + A2'/(2 Lam_{i+1} s^3) - Lam_{i+1}^2 N(beta; Lam_i)
     q_i = Lam_{i+1}/Lam_i is chosen as large as allowed by  r < 1 on beta >= 0.37  and  rt < 1 on beta >= 0.245  (cells + tail).
  Segment B (Lam >= Lam_a, anchors Lam_a): with the envelope R^ (gl_common), Lam^4 N^(beta; Lam) is non-increasing [P], so
        r <= (Lam_a/Lam)^2 r^(beta; Lam_a),   rt <= (Lam_a/Lam) rt^(beta; Lam_a),   Lam^2 F_LR >= min(A1min M1, Lam_a^2 F^_cert(Lam_a)).
     Anchor Lam_n closes the range; anchors Lam(5000), Lam(1e4), Lam(2e4), Lam(5e4) give the Dg-improved thresholds beta0(Dg).
Everything in Arb (128 bits), decisions by Arb comparisons.  Output: cache/gl2_allgaps.{json,txt}
"""
import os, sys, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from flint import arb
import gl_common as G
from gl_common import A, up, lo, W, A1MIN, A2P

lines = []
def pr(s=''):
    print(s, flush=True); lines.append(s)

B_R = Fr(37, 100); B_T = Fr(245, 1000)
QCAND = [Fr(n, 10000) for n in (10005, 10010, 10015, 10020, 10030, 10040, 10050, 10075, 10100, 10150, 10200, 10300, 10400, 10500,
                                10750, 11000, 11500, 12000, 13000, 14000, 15000)]


def cell_data(LAM, cl, hat):
    TAB = {k: G.build_I(W[k], LAM, hat=hat) for k in (2, 3)}
    out = []
    for (a, b, m1) in cl:
        e2, _ = G.E_up(2, a, b, TAB, LAM); e3, _ = G.E_up(3, a, b, TAB, LAM)
        out.append(dict(a=a, b=b, A1M1=lo(A(A1MIN)*m1), L2N=up(LAM**2*(2*e2 + e3)), sb3=up((2*(A(b)/2).sinh())**3)))
    tl, _ = G.tail_majorant(TAB, LAM)        # = Lam^2 N/(A1min M1) bound for beta >= 64
    return out, tl


def interval_bounds(cd, tl, q, Lb):
    """bounds valid for Lam in [Lam_i, q Lam_i], from data cd computed at Lam_i (L2N = Lam_i^2 N(Lam_i)); Lb = q Lam_i."""
    q2 = A(q)**2
    res = []
    for c in cd:
        num = q2*c['L2N']
        flat = A(A2P)/(2*A(Lb)*c['sb3'])                 # A2'/(2 Lam s^3) is decreasing in Lam: value at the right end
        res.append(dict(a=c['a'], r=up(num/c['A1M1']), rt=up(num/(c['A1M1'] + flat)), L2F=lo(c['A1M1'] + flat - num)))
    return res, up(q2*tl)


def ok(res, tlq):
    if not tlq < 1:
        return False
    for x in res:
        if x['a'] >= B_R and not x['r'] < 1:
            return False
        if x['a'] >= B_T and not x['rt'] < 1:
            return False
    return True


def sup_key(res, key, bp):
    v = arb(0)
    for x in res:
        if x['a'] >= bp and x[key].upper() > v.upper():
            v = x[key]
    return v


def main():
    pr(__doc__.split('Output')[0])
    t0 = time.time()
    cl = G.m1_cells()
    floor_min = {}                                       # per cell: min over all Lam >= Lam0 of the gap-unit floor bound
    def upd(res):
        for x in res:
            k = x['a']
            if k not in floor_min or x['L2F'].upper() < floor_min[k].upper():
                floor_min[k] = x['L2F']
    # ---------------- segment A
    Lam_end = G.lam_of(2000)
    Li = Fr(998499874812, 10**9) - Fr(1, 10**6)          # rational lower bound of Lam(1000) = 998.499874812...
    assert A(Li) < G.LAM0
    grid = []
    pr('Segment A (LRINF_TIGHT bound with R, Lam non-increasing pieces):')
    pr('  i | Lam_i | q_i | sup r (beta >= 0.37) | sup rt (beta >= 0.245) on [Lam_i, Lam_i+1] | min over cells beta >= 0.37 of (Lam^2 F bound)/(A1min M1)')
    while not A(Li) > Lam_end:
        cd, tl = cell_data(A(Li), cl, hat=False)
        best = None
        for q in QCAND:
            res, tlq = interval_bounds(cd, tl, q, Li*q)
            if ok(res, tlq):
                best = (q, res, tlq)
            else:
                break
        assert best is not None, f'no admissible step at Lam = {float(Li)}'
        q, res, tlq = best
        upd(res)
        sr = sup_key(res, 'r', B_R); st = sup_key(res, 'rt', B_T)
        A1 = {c['a']: c['A1M1'] for c in cd}
        mf = min((x['L2F']/A1[x['a']] for x in res if x['a'] >= B_R), key=lambda v: float(v.mid()))
        grid.append(dict(Lam=float(Li), q=float(q), sup_r=float(sr.mid()), sup_rt=float(st.mid()), tail=float(tlq.mid())))
        pr(f'  {len(grid)-1:2d} | {float(Li):10.3f} | {float(q):.4f} | {float(sr.mid()):.4f} | {float(st.mid()):.4f} | {float(mf.mid()):.4f}'
           f'   ({time.time()-t0:.0f}s)')
        Li = Li*q
    LamN = Li
    pr(f'  segment A covers [{grid[0]["Lam"]:.3f}, {float(LamN):.3f}] in {len(grid)} steps (Lam(2000) = {Lam_end.str(10)})')
    # ---------------- segment B anchors
    out = dict(segA=grid, LamN=float(LamN))
    anchors = [('Lam_n', A(LamN), float(LamN))] + [(f'Lam({Dg})', G.lam_of(Dg), Dg) for Dg in (5000, 10000, 20000, 50000)]
    pr('\nSegment B anchors (envelope R^, Lam^4 N^ non-increasing): bounds valid for every Lam >= Lam_a')
    anc = []
    for name, La, tag in anchors:
        cd, tl = cell_data(La, cl, hat=True)
        res = []
        for c in cd:
            num = c['L2N']; flat = A(A2P)/(2*La*c['sb3'])
            res.append(dict(a=c['a'], r=up(num/c['A1M1']), rt=up(num/(c['A1M1'] + flat)),
                            L2F=c['A1M1'] if c['A1M1'] < c['A1M1'] + flat - num else lo(c['A1M1'] + flat - num)))
        assert tl < 1
        def b0(key):
            v = None
            for x in reversed(res):
                if x[key] < 1:
                    v = x['a']
                else:
                    break
            return v
        b0r, b0t = b0('r'), b0('rt')
        sr = sup_key(res, 'r', B_R); st = sup_key(res, 'rt', B_T)
        d = dict(name=name, Lam=float(La.mid()), beta0_r=str(b0r), beta0_rt=str(b0t), sup_r_037=float(sr.mid()),
                 sup_rt_0245=float(st.mid()), tail=float(tl.mid()))
        anc.append(d)
        pr(f'  {name:12s} Lam_a = {float(La.mid()):10.3f}: sup r^ (beta >= 0.37) = {d["sup_r_037"]:.4f}; sup rt^ (beta >= 0.245) = '
           f'{d["sup_rt_0245"]:.4f}; beta0(r^ < 1) = {float(b0r):.4f}; beta0(rt^ < 1) = {float(b0t):.4f}; tail {d["tail"]:.2e}')
        if name == 'Lam_n':
            assert sr < 1 and st < 1
            upd(res)
            closing = res
    out['anchors'] = anc
    # ---------------- uniform floor table
    pr('\nuniform gap-unit floor over ALL Dg >= 1000:  Lam^2 F_LR(beta, u, sigma) >= f_U(cell)  (u <= 1, sigma = +-1)')
    fl = []
    for c in cl:
        a = c[0]
        v = floor_min[a]
        fl.append(dict(a=str(a), b=str(c[1]), fU=float(v.mid()), A1M1=float(lo(A(A1MIN)*c[2]).mid())))
    for a in ('0.245', '0.25', '0.3', '0.35', '0.37', '0.4', '0.5', '0.7', '1', '1.5', '2', '3', '5', '10', '20', '40'):
        x = next(f for f in fl if Fr(f['a']) == Fr(a))
        pr(f"  [{float(Fr(x['a'])):.3f}, {float(Fr(x['b'])):.3f}]: f_U >= {x['fU']:.5f}   (A1min M1 >= {x['A1M1']:.5f}; ratio {x['fU']/x['A1M1']:.3f})")
    fpos = all(Fr(f['a']) < B_T or floor_min[Fr(f['a'])] > 0 for f in fl)
    pr(f'  f_U > 0 on every cell with beta >= 0.245: {fpos}')
    for thr in (Fr(99, 100), Fr(9, 10), Fr(1, 2)):
        bf = None
        for c in reversed(cl):
            if floor_min[c[0]] >= A(thr)*lo(A(A1MIN)*c[2]):
                bf = c[0]
            else:
                break
        pr(f'  smallest beta_f with f_U >= {float(thr)} x A1min M1 on every cell from beta_f on: {float(bf)}')
        out[f'beta_f_{float(thr)}'] = str(bf)
    out['floor'] = fl
    json.dump(out, open(os.path.join(G.CACHE, 'gl2_allgaps.json'), 'w'), indent=0)
    open(os.path.join(G.CACHE, 'gl2_allgaps.txt'), 'w').write('\n'.join(lines) + '\n')
    pr(f'done in {time.time()-t0:.0f}s')


if __name__ == '__main__':
    main()
