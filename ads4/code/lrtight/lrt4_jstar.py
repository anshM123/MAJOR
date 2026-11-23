"""
lrtight/lrt4_jstar.py (LRINF_TIGHT.md §3): J*(m) for a Lemma LR-inf threshold beta0, decided exactly / in Arb.

  beta(J, m) >= beta0  <=>  bhat = (2J+1) Lam/m >= bhat0 = 2 Lam sinh(beta0/2)  <=>  (2J+1) >= s0 m,  s0 = 2 sinh(beta0/2).
  * fixed m: J*(m) = ceil((s0 m - 1)/2);
  * twist edge m^2 = m_gap(J)^2 = (tau-2)(tau+2J-1) (CHMRS d = 3, tau = 1000), and c * m_gap: (2J+1)^2 >= s0^2 c^2 (tau-2)(tau+2J-1);
    (2J+1)^2/(tau+2J-1) is increasing in J, so J* is the first J satisfying it and every J >= J* does.
  The comparisons are made in Arb; an undecidable comparison raises.
Rows: beta0 = 1 (reproduces LARGESPIN_L: 857, 987, 1275, 780, 1041, 2602), the certified thresholds of lrt3 (r < 1, rtilde < 1,
rtilde < 1/2), and beta0 = 0.5 (INDEP_CHECK2's side result, 324).
Output: cache/lrt4_jstar.json, cache/lrt4_jstar.txt
"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from flint import arb, ctx
from lrt_common import A, LAM, CACHE, TAU

ctx.prec = 128
lines = []
def pr(s=''):
    print(s, flush=True); lines.append(s)


def geq(x, y):
    d = x - y
    if d >= 0:
        return True
    if d < 0:
        return False
    raise ValueError('undecided comparison')


def jstar_gap(beta0, c):
    s0 = 2*(A(beta0)/2).sinh()
    c = A(c)
    J = 0
    while not geq(A((2*J + 1)**2), s0**2*c**2*(TAU - 2)*(TAU + 2*J - 1)):
        J += 1
    return J


def jstar_fixed(beta0, mL):
    s0 = 2*(A(beta0)/2).sinh()
    v = (s0*A(mL)*LAM - 1)/2
    f = v.floor(); cl = v.ceil()
    assert f.is_exact() and cl.is_exact(), 'ceil undecided'
    return int(cl.unique_fmpz())


ROWS_M = [('m_gap', 'g', Fr(1)), ('1.1 m_gap', 'g', Fr(11, 10)), ('1.3 m_gap', 'g', Fr(13, 10)),
          ('1.5 Lam', 'f', Fr(3, 2)), ('2 Lam', 'f', Fr(2)), ('5 Lam', 'f', Fr(5))]


def main():
    pr(__doc__.split('Output')[0])
    R = json.load(open(os.path.join(CACHE, 'lrt3_rcert.json')))
    cases = [('beta0 = 1 (LARGESPIN_L, r <= 0.577 crude / 0.0162 tight)', Fr(1)),
             ('beta0 = 0.5 (INDEP_CHECK2 side result)', Fr(1, 2)),
             ('beta0(r < 1/2), certified', Fr(R['beta0_r_0.5'])),
             ('beta0(r < 1), certified  [Lemma LR-inf as stated]', Fr(R['beta0_r_1.0'])),
             ('beta0(rtilde < 1/2), certified', Fr(R['beta0_rt_0.5'])),
             ('beta0(rtilde < 1), certified  [F_LR > 0]', Fr(R['beta0_rt_1.0']))]
    out = {}
    pr('| threshold | beta0 | bhat0/Lam | ' + ' | '.join(r[0] for r in ROWS_M) + ' |')
    for name, b0 in cases:
        s0 = 2*(A(b0)/2).sinh()
        js = []
        for lab, kind, c in ROWS_M:
            js.append(jstar_gap(b0, c) if kind == 'g' else jstar_fixed(b0, c))
        out[str(b0)] = dict(name=name, s0=float(s0.mid()), J=js)
        pr(f'| {name} | {float(b0):.3f} | {float(s0.mid()):.5f} | ' + ' | '.join(str(j) for j in js) + ' |')
    # m_gap at the new J* (which m-band the gap states fall in)
    pr('\nm_gap(J)/Lam at the twist-edge J*: ' + ', '.join(
        f"{out[k]['J'][0]}: {float((A((TAU - 2)*(TAU + 2*out[k]['J'][0] - 1)).sqrt()/LAM).mid()):.4f}" for k in out))
    pr('formula: J*(m) = ceil((s0 m - 1)/2), s0 = 2 sinh(beta0/2) (m in AdS units).')
    json.dump(out, open(os.path.join(CACHE, 'lrt4_jstar.json'), 'w'), indent=1)
    open(os.path.join(CACHE, 'lrt4_jstar.txt'), 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
