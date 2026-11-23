"""
lrtight/lrt3_rcert.py (LRINF_TIGHT.md §2.3-2.5): certified upper bound of Lemma LR-inf's remainder ratio on beta-cells covering
[0.1, 64], plus a certified tail bound for beta >= 64.

For a cell [ba, bb] and k = 2, 3 (every piece below is monotone in beta, proved in LRINF_TIGHT §2.3):
   Ta  <= K0(ba)/Lam * d1(Lam ba)                       (K0 = (cosh b - 1)^{-1/2} and d1 decreasing)
   Tb  <= K0(ba)/(4 D(ba) Lam^3) * sumc(Lam ba)          (D = cosh b - 1 increasing, sumc decreasing)
   Tc1 <= K0(ba) [ I1(bb/2)/(48 D(ba)) + (3/8) 2^{5/2} I2(bb/2)/D(ba)^2 ],
          I1(T) = int_0^T R(Lam t) t^4 cosh t dt,  I2(T) = int_0^T R(Lam t)(cosh t - 1)^2 dt   (rho_b decreasing in beta; range grows)
   Tc2 <= R_inf(Lam ba/2) (pi/sqrt2) P_{-1/2}(cosh ba)   (R_inf, P_{-1/2}(cosh b) = 2F1(1/2,1/2;1;-sinh^2(b/2)) decreasing)
   |E_k| <= (sqrt2/pi)(Ta + Tb + Tc1 + Tc2)
I1, I2: rigorous upper Riemann sums on a geometric t-grid (R(Lam t) decreasing, t^4 cosh t and (cosh t - 1)^2 increasing), plus the
explicit [0, t1] piece  K: int_0^{t1} (a1/(Lam t)^2 + int|psi|) t^4 C dt.
M1: certified lower bounds of lrt2_m1 (cells, beta <= 24) and the analytic region-(c) bound evaluated on the cell ball (beta > 24).
Outputs per cell (u <= 1, both sigma, i.e. all m >= Lam and both helicities):
   r      = (2|E2| + |E3|)/(A1min M1/Lam^2)                         (LARGESPIN_L's r(beta), tight constant int|psi|)
   rtilde = (2|E2| + |E3|)/(A1min M1/Lam^2 + A2'/(2 bhat^3))         (the flat tail included; F_LR > 0 iff-sufficient: rtilde < 1)
   Fcert  = A1min M1/Lam^2 + A2'/(2 bhat^3) - (2|E2| + |E3|)         (certified lower bound of min_{u<=1,sigma} F_LR on the cell)
Tail beta >= 64: explicit majorant, a sum of decreasing functions (LRINF_TIGHT §2.4), evaluated at 64 in Arb.
Output: cache/lrt3_rcert.json, cache/lrt3_rcert.txt
"""
import os, sys, json, time, math
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from flint import arb, ctx
from lrt_common import Weight, A, LAM, CACHE, up, lo, ball, P_mhalf
from lrt2_m1 import M1_lower_region_c, cells

ctx.prec = 128
lines = []
def pr(s=''):
    print(s, flush=True); lines.append(s)

CONST = json.load(open(os.path.join(CACHE, 'lrt1_constants.json')))
W = {k: Weight(k, absint_up=Fr(CONST[str(k)]['absint_upper_rational'])) for k in (2, 3)}
A1MIN = Fr(CONST['A1min']); A2P = Fr(CONST['A2p'])
C3 = arb(3)/8*arb(2)**(arb(5)/2)
SQ2PI = arb(2).sqrt()/arb.pi()

# ---------------- I1, I2 cumulative upper sums ----------------
T1 = 1e-6; Q = 1.0005; TMAX = 32.5


def build_I(w):
    ts = [T1]
    while ts[-1] < TMAX:
        ts.append(ts[-1]*Q)
    ts = [arb(t) for t in ts]                     # exact dyadic partition points
    t1 = ts[0]
    CC = t1.cosh()
    a1 = A(w.a1); ai = A(w.absint_up)
    I1 = a1*t1**3/(3*LAM**2)*CC + ai*t1**5/5*CC            # [0, t1]: R <= a1/x^2 + int|psi|, t^4 cosh t <= t^4 cosh t1
    I2 = (a1*t1**3/(3*LAM**2) + ai*t1**5/5)*CC**2/4         # (cosh t - 1)^2 <= t^4 cosh^2(t1)/4
    cum1 = [up(I1)]; cum2 = [up(I2)]
    for i in range(len(ts) - 1):
        a, b = ts[i], ts[i + 1]
        Ra = w.R(LAM*a)
        h = b - a
        I1 = I1 + up(Ra*b**4*b.cosh()*h)
        I2 = I2 + up(Ra*(b.cosh() - 1)**2*h)
        cum1.append(up(I1)); cum2.append(up(I2))
    return ts, cum1, cum2


def lookup(ts, cum, T):
    """upper bound of int_0^T (nonneg integrand): first grid point >= T."""
    lo_, hi_ = 0, len(ts) - 1
    assert ts[-1] >= T
    while lo_ < hi_:
        mid = (lo_ + hi_)//2
        if ts[mid] >= T:
            hi_ = mid
        else:
            lo_ = mid + 1
    return cum[lo_]


def E_up(k, ba, bb, TAB):
    w = W[k]; ts, c1, c2 = TAB[k]
    Ba = A(ba); Bb = A(bb)
    D = Ba.cosh() - 1; K0 = 1/D.sqrt(); X = LAM*Ba
    Ta = K0/LAM*w.d1f(X)
    Tb = K0/(4*D*LAM**3)*w.sumc(X)
    Tc1 = K0*(lookup(ts, c1, Bb/2)/(48*D) + C3*lookup(ts, c2, Bb/2)/D**2)
    Tc2 = w.Rinf(LAM*Ba/2)*(arb.pi()/arb(2).sqrt())*P_mhalf(Ba)
    E = SQ2PI*(Ta + Tb + Tc1 + Tc2)
    return up(E), dict(Ta=SQ2PI*Ta, Tb=SQ2PI*Tb, Tc1=SQ2PI*Tc1, Tc2=SQ2PI*Tc2)


def tail_majorant(TAB, BT=64):
    """r(beta) <= value of a sum of decreasing functions at BT, for every beta >= BT (LRINF_TIGHT §2.4)."""
    B = arb(BT); T0 = arb(32)
    ln2 = arb(2).log()
    ell = 1 - B/(2*(B/2).sinh()) - arb(2).sqrt()*(B*(-B/2).exp()/A(Fr(47, 10)) + 4/B) - 2*(2*ln2).sqrt()*B/(B - ln2)**2   # beta L(beta), increasing
    assert ell > 0
    D = B.cosh() - 1
    C1 = 2*arb(2).sqrt()*(arb(1)/4).exp(); C2 = (2/(1 - (-arb(1)).exp())).sqrt()
    tot = arb(0); parts = {}
    for k, wk in ((2, 2), (3, 1)):
        w = W[k]; ts, c1, c2 = TAB[k]
        I1T0 = lookup(ts, c1, T0); I2T0 = lookup(ts, c2, T0)
        rinf0 = w.Rinf(LAM*T0)
        rho0 = rinf0*T0**4
        f_a = B*w.d1f(LAM*B)/LAM
        f_b = B*w.sumc(LAM*B)/(4*D*LAM**3)
        f_c1 = B*I1T0/(48*D) + rho0*B/(96*(B/2).sinh())
        f_c2 = C3*(B*I2T0/D**2 + rinf0*(B*B/(16*(B/2).sinh()**4) + B*(B/2).cosh()/(8*(B/2).sinh()**3)))
        f_d = B*w.Rinf(LAM*B/2)*(C1 + C2*B)/arb(2).sqrt()
        s = f_a + f_b + f_c1 + f_c2 + f_d
        parts[k] = s
        tot += wk*s
    r = LAM**2/(A(A1MIN)*ell)*tot
    return up(r), parts, ell


def main():
    pr(__doc__.split('Output')[0])
    t0 = time.time()
    TAB = {k: build_I(W[k]) for k in (2, 3)}
    pr(f'I1, I2 grids built ({len(TAB[2][0])} points, {time.time()-t0:.0f}s): I1(1/2) k=2 <= {TAB[2][1][len(TAB[2][0])//2].str(5)}')
    M1c = json.load(open(os.path.join(CACHE, 'lrt2_m1.json')))
    cl = [(Fr(a), Fr(b), arb(m)) for a, b, m in M1c]
    # beta in (24, 64]: cells of 0.5 with the region (c) bound on the cell ball
    for (a, b) in cells([('24', '64', '0.5')]):
        cl.append((a, b, lo(M1_lower_region_c(ball(a, b)))))
    rows = []
    for (a, b, m1) in cl:
        assert m1 > 0
        e2, p2 = E_up(2, a, b, TAB); e3, p3 = E_up(3, a, b, TAB)
        num = 2*e2 + e3
        main_ = A(A1MIN)*m1/LAM**2
        bh_b = 2*LAM*(A(b)/2).sinh()
        flat = A(A2P)/(2*bh_b**3)
        r = up(num/main_); rt = up(num/(main_ + flat)); Fc = lo(main_ + flat - num)
        rows.append(dict(ba=str(a), bb=str(b), M1lo=float(m1.mid()), E2=float(e2.mid()), E3=float(e3.mid()),
                         r=float(r.mid()), rt=float(rt.mid()), Fcert=float(Fc.mid()),
                         parts2={k: float(up(v).mid()) for k, v in p2.items()}))
    rT, partsT, ell = tail_majorant(TAB)
    pr(f'tail: for every beta >= 64, r(beta) <= {rT.str(5)}  (beta L(beta) >= {ell.str(6)} at 64; per-k majorants {partsT[2].str(4)}, {partsT[3].str(4)})')
    # summary
    pr('\ncell table (selected): beta-cell | M1 lower | 2|E2|+|E3| upper | r upper | rtilde upper | F_cert (min over u<=1, sigma)')
    for R_ in rows:
        a = Fr(R_['ba'])
        if (a*20).denominator == 1 and (a <= 3 or a.denominator == 1):
            pr(f"  [{float(a):6.3f}, {float(Fr(R_['bb'])):6.3f}]  M1 >= {R_['M1lo']:.5e}  2E2+E3 <= {2*R_['E2']+R_['E3']:.3e}  "
               f"r <= {R_['r']:.4f}  rtilde <= {R_['rt']:.4f}  Fcert = {R_['Fcert']:+.4e}   [k=2: Ta {R_['parts2']['Ta']:.1e} Tb {R_['parts2']['Tb']:.1e} Tc1 {R_['parts2']['Tc1']:.1e} Tc2 {R_['parts2']['Tc2']:.1e}]")
    # sup over beta >= beta'
    def supfrom(key, bp):
        return max([R_[key] for R_ in rows if Fr(R_['ba']) >= bp - Fr(1, 10**9)] + [float(rT.mid())])
    pr('\nsup over beta >= beta\' (cells to 64 + tail):')
    for bp in ['0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', '0.45', '0.5', '0.6', '0.7', '0.8', '0.9', '1', '1.2', '1.5', '2', '3', '5', '24', '64']:
        bp = Fr(bp)
        if bp == 64:
            pr(f"  beta >= 64: r <= {float(rT.mid()):.3e} (tail)")
            continue
        pr(f"  beta >= {float(bp):5.2f}: sup r <= {supfrom('r', bp):.4e}   sup rtilde <= {supfrom('rt', bp):.4e}")
    # beta0: smallest cell start such that every cell from there on has r < 1 (resp. rtilde < 1)
    def beta0(key, thr=1.0):
        b0 = None
        for R_ in reversed(rows):
            if R_[key] < thr:
                b0 = Fr(R_['ba'])
            else:
                break
        return b0
    out = dict(rows=rows, tail=float(rT.mid()))
    for key in ('r', 'rt'):
        for thr in (1.0, 0.5):
            b0 = beta0(key, thr)
            out[f'beta0_{key}_{thr}'] = str(b0)
            pr(f'smallest certified beta0 with {key} < {thr} on [beta0, oo): {float(b0) if b0 is not None else None}')
    json.dump(out, open(os.path.join(CACHE, 'lrt3_rcert.json'), 'w'), indent=0)
    open(os.path.join(CACHE, 'lrt3_rcert.txt'), 'w').write('\n'.join(lines) + '\n')
    pr(f'done in {time.time()-t0:.0f}s')


if __name__ == '__main__':
    main()
