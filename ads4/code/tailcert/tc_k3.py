"""
tailcert/tc_k3.py (TAIL_CERT.md §5): the k = 3 piece H^(3) = w^3 sigma_b tau L(G_D), G_D = Lam^2 int_0^1 psi3(p) F_{4,Lam p} dp
(FULLPOS3 §4.1(a); no dressing on G_D), bounded in absolute value on the real frequency axis:
  |L(F_{4,nu})| <= 2|B(nu)| w^{-7/2} [ m0 + sum_{j>=1} |c_j(nu)| w^j m_j ] + |L(R'_{J0})|,
  m0 >= |Ŝ(7/2 + i nu)|, m_{-1} >= |Ŝ(9/2 + i nu)| (tc_power Theorem P), m_j by the exact recurrence (tc_exact E3).
Relative to the main term (2/Lam) sqrt(w) I_W(L):
  |H^(3)|/M <= (tau/w) Lam^3 int_0^1 |psi3||B(Lam p)| X3(p) dp / I_W(L) + (remainder),   tau/w = 1/(T-C) <= 0.2515 w/a^2.
Usage: python3 tc_k3.py L0 J0
"""
import os, sys, json, math
from fractions import Fraction as Fr
from flint import arb, acb, ctx
ctx.prec = 128
HERE = os.path.dirname(os.path.abspath(__file__))
ADS4 = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)
import tc_power as P
import tc_assemble as AS
import tc_remainder as RM
LAM = arb(998*999).sqrt()


def up(x):
    return float(arb(x).upper())


def main():
    L0 = float(sys.argv[1]); J0 = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    w = math.exp(-L0)
    out = []
    def pr(s):
        print(s, flush=True); out.append(s)
    nf = json.load(open(os.path.join(ADS4, 'prl', 'data', 'nf02_psi.json')))
    K0c = Fr('27.48509817')
    p3 = [Fr(c)/K0c for c in nf['psi3_exact']]
    assert p3[0] == 0
    import math as _m
    from flint import fmpq_poly, fmpq
    q3poly = fmpq_poly([fmpq(c.numerator, c.denominator) for c in p3[1:]])     # psi3(p)/p, exact
    NC = 2000
    tot = arb(0); totabs = arb(0)
    for k in range(NC):
        p1, p2_ = k/NC, (k + 1)/NC
        pb = arb((p1 + p2_)/2, (p2_ - p1)/2)
        # sup |psi3/p| on the cell: exact Taylor shift to the cell midpoint (fmpq), then sum |d_k| r^k
        pm = fmpq(2*k + 1, 2*NC); rr = arb(1)/(2*NC)
        sh = q3poly(fmpq_poly([pm, 1]))
        psi_over_p = sum((abs(arb(sh[i].p)/arb(sh[i].q))*rr**i for i in range(sh.degree() + 1)), arb(0))
        nu = acb(LAM*pb)
        # |nu B(nu)|^2 = nu coth(pi nu)/(4 pi) for real nu (increasing): sup at the right end of the cell
        n2 = LAM*p2_
        nuB = arb(((n2*(arb.pi()*n2).coth())/(4*arb.pi())).sqrt().upper())
        t1, t2 = float((LAM*p1).lower()), float((LAM*p2_).upper())
        m0, _ = P.K_osc(3.5, max(t1, 0.0), t2, L0)
        mm1, _ = P.K_osc(4.5, max(t1, 0.0), t2, L0)
        cbs = AS.cj_bounds(nu, J0)
        # recurrence without dressing: reuse jsum with w -> dressing negligible (it only adds); m_j bound identical
        js, X0 = AS.jsum(m0, mm1, abs(complex(3.5, t2)), cbs, w, J0)
        X3 = X0 + js
        tot += psi_over_p*nuB/LAM*arb(X3)*(arb(p2_) - p1)
    pr(f'k=3: int_0^1 |psi3(p)| |B(Lam p)| X3(p) dp <= {up(tot):.4e}   (L0 = {L0}, J0 = {J0})')
    W_I = json.load(open(os.path.join(ADS4, 'code', 'contour', 'cache', 'cb_w.json')))['cells']
    worst = 0.0
    for L in [L0 + 0.5*k for k in range(200)]:
        IW = sum(mn*(math.exp(-L*y0) - math.exp(-L*y1))/L for y0, y1, mn, mx_ in W_I)
        rel = 0.2515*math.exp(-L)*up(LAM**3)/998**2*up(tot)/IW
        worst = max(worst, rel)
    pr(f'k=3 main part / M(L) <= {worst:.3e} for L in [{L0}, {L0 + 99.5}] (sup over the grid; the ratio is decreasing in L for L > 1)')
    # remainder (p = 4, e0 = 7/2, m = J0 - 5): |L(R')| <= 2.53 K (1.0006 w)^{m-1} Q_m (1 + 40*1.5^m)  (dressing-free version is smaller)
    m = J0 - 5
    lK0, cL = RM.logK_small(4, J0, m, L0)
    lsum = []
    for k in range(1, 1000):
        n1, n2 = float((LAM*k/1000).lower()), float((LAM*(k + 1)/1000).upper())
        if n1 < 0.4:
            continue
        lK = RM.logK_real(n1, n2, 4, J0, m)
        lsum.append(lK + math.log(1/1000))
    mx = max(lsum); lI = mx + math.log(sum(math.exp(x - mx) for x in lsum))
    # sup|psi3| <= sum |coeff|
    s3 = sum(abs(float(c)) for c in p3)
    worst_r = -1e300
    for L in [L0 + 0.5*k for k in range(200)]:
        IW = sum(mn*(math.exp(-L*y0) - math.exp(-L*y1))/L for y0, y1, mn, mx_ in W_I)
        lrel = math.log(0.2515*s3) - L + 3*math.log(float(LAM.mid())) - 2*math.log(998) + 2.5*(-L) + RM.log_loc_factor(m, L) \
            + max(lI, lK0 + cL*L) + math.log(2) - math.log(IW) + math.log(float(LAM.mid()))
        worst_r = max(worst_r, lrel)
    pr(f'k=3 remainder / M(L): log10 <= {worst_r/math.log(10):.1f}')
    open(os.path.join(HERE, 'cache', f'tc_k3_L{L0}_J{J0}.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
