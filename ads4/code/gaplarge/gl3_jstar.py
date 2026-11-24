"""
gaplarge/gl3_jstar.py (GAP_LARGESPIN.md §1.6): J*(m; Dg) for the Dg-uniform and the Dg-improved certified thresholds of gl2, decided in
Arb exactly as lrtight/lrt4_jstar.py:
  beta >= beta0  <=>  2J + 1 >= s0 m  (m in AdS units),  s0 = 2 sinh(beta0/2);   fixed m = mu Lam: J* = ceil((s0 mu Lam - 1)/2);
  twist edge m^2 = (Dg-2)(Dg+2J-1) (and c x m_gap): first J with (2J+1)^2 >= s0^2 c^2 (Dg-2)(Dg+2J-1) (LHS/(Dg+2J-1) increasing in J).
  Asymptotics: J*_edge/Lam -> x*(s0) = (s0^2 + s0 sqrt(s0^2 + 4))/4 as Dg -> oo  (from 4x^2 = s0^2 (1 + 2x)).
Output: cache/gl3_jstar.{txt,json}
"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from flint import arb, ctx
ctx.prec = 128
import gl_common as G
from gl_common import A

lines = []
def pr(s=''):
    print(s, flush=True); lines.append(s)


def geq(x, y):
    d = x - y
    if d >= 0:
        return True
    if d < 0:
        return False
    raise ValueError('undecided')


def jstar_gap(b0, c, Dg):
    s0 = 2*(A(b0)/2).sinh(); c = A(c)
    # start near the asymptotic solution and step down/up exactly
    J = 0
    lo_, hi_ = 0, 10**9
    while lo_ < hi_:                 # monotone predicate: bisection
        mid = (lo_ + hi_)//2
        if geq(A((2*mid + 1)**2), s0**2*c**2*(Dg - 2)*(Dg + 2*mid - 1)):
            hi_ = mid
        else:
            lo_ = mid + 1
    return lo_


def jstar_fixed(b0, mu, Dg):
    s0 = 2*(A(b0)/2).sinh()
    v = (s0*A(mu)*G.lam_of(Dg) - 1)/2
    cl = v.ceil()
    assert cl.is_exact()
    return int(cl.unique_fmpz())


COLS = [('m_gap', 'g', Fr(1)), ('1.1 m_gap', 'g', Fr(11, 10)), ('1.3 m_gap', 'g', Fr(13, 10)),
        ('1.5 Lam', 'f', Fr(3, 2)), ('2 Lam', 'f', Fr(2)), ('5 Lam', 'f', Fr(5))]


def main():
    pr(__doc__.split('Output')[0])
    g2 = json.load(open(os.path.join(HERE, 'cache', 'gl2_allgaps.json')))
    an = {a['name']: a for a in g2['anchors']}
    LamN = g2['LamN']; DgN = 1.5 + (LamN**2 + 0.25)**0.5
    pr(f'segment-A/anchor boundary: Lam_n = {LamN:.3f}  <=>  Dg_n = {DgN:.2f}')
    regimes = [(1000, 'all Dg >= 1000 (uniform)', Fr(37, 100), Fr(245, 1000)),
               (2606, f'Dg >= 2606 (anchor Lam_n)', Fr(an['Lam_n']['beta0_r']), Fr(an['Lam_n']['beta0_rt'])),
               (5000, 'Dg >= 5000', Fr(an['Lam(5000)']['beta0_r']), Fr(an['Lam(5000)']['beta0_rt'])),
               (10000, 'Dg >= 1e4', Fr(an['Lam(10000)']['beta0_r']), Fr(an['Lam(10000)']['beta0_rt'])),
               (20000, 'Dg >= 2e4', Fr(an['Lam(20000)']['beta0_r']), Fr(an['Lam(20000)']['beta0_rt']))]
    out = {}
    pr('\n(a) Dg-uniform thresholds evaluated at several Dg (r < 1 form: beta0 = 0.37, s0 = 0.37211; F_LR > 0: beta0 = 0.245):')
    pr('| Dg | Lam | form | ' + ' | '.join(c[0] for c in COLS) + ' | J*_edge/Lam |')
    for Dg in (1000, 1500, 2000, 3000, 5000, 10000, 100000):
        L = float(G.lam_of(Dg).mid())
        for form, b0 in (('r<1', Fr(37, 100)), ('F>0', Fr(245, 1000))):
            js = [jstar_gap(b0, c, Dg) if k == 'g' else jstar_fixed(b0, c, Dg) for _, k, c in COLS]
            out[f'{Dg}|{form}|uniform'] = js
            pr(f'| {Dg} | {L:.1f} | {form} | ' + ' | '.join(str(j) for j in js) + f' | {js[0]/L:.4f} |')
    for b0 in (Fr(37, 100), Fr(245, 1000)):
        s0 = float((2*(A(b0)/2).sinh()).mid())
        pr(f'   Dg -> oo: J*_edge/Lam -> x*(s0 = {s0:.5f}) = {(s0*s0 + s0*(s0*s0 + 4)**0.5)/4:.5f};  J*(mu Lam)/Lam -> {s0/2:.5f} mu')
    pr('\n(b) Dg-improved certified thresholds (valid for every Dg in the stated half-line; beta-cells stop at 0.1):')
    pr('| regime | beta0 (r<1) | s0 | beta0 (F>0) | s0 | J* at Dg = regime start, r<1 form: ' + ' | '.join(c[0] for c in COLS) + ' |')
    for Dg, lab, br, bt in regimes:
        sr = float((2*(A(br)/2).sinh()).mid()); st = float((2*(A(bt)/2).sinh()).mid())
        js = [jstar_gap(br, c, Dg) if k == 'g' else jstar_fixed(br, c, Dg) for _, k, c in COLS]
        out[f'{Dg}|improved'] = dict(beta0_r=str(br), beta0_rt=str(bt), J=js)
        pr(f'| {lab} | {float(br):.4f} | {sr:.5f} | {float(bt):.4f} | {st:.5f} | ' + ' | '.join(str(j) for j in js) + ' |')
    json.dump(out, open(os.path.join(HERE, 'cache', 'gl3_jstar.json'), 'w'), indent=0)
    open(os.path.join(HERE, 'cache', 'gl3_jstar.txt'), 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
