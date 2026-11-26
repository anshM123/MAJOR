"""
rt_loc3.py: Lemma Loc' (ANALYTIC_SLIVER §7.3) -- localised bound for the j-truncation remainder, replacing TAIL_CERT's Lemma Loc in
Route T' below L = 8.  Supersedes the probe rt_loc2.py (adds: dressing, tiered B-regions, small-nu pole handling, nu-integration).

Setting (TAIL_CERT §4).  R = R_J(xi; nu) = F_{p,nu}(xi) - T_J,  T_J = sum_{+-} B(+-nu) sum_{j<J} c_j(+-nu) (1-xi)^{j-e0 -+ i nu},  e0 = p - 1/2,
analytic in |xi| < 1; on |1 - xi| < 1 also R = sum_{+-} B(+-nu) sum_{j>=J} c_j(+-nu)(1-xi)^{j-e0-+i nu}.
  L(R) = E_B[ (1/2 pi i) oint_{|z|=rho} (1-z)^{-a-1} R(Bz) z^{-n-1} dz ] / C(N,a),   B ~ Beta(n+s, a),  rho = n/N.
  K0 = rho^{-n} I_W / C(N,a) <= 1.6,  I_W = (1/2pi) int |1 - rho e^{i th}|^{-a-1} dth  (Stirling/Robbins, TAIL's I_W bound, Gautschi);
  fraction of the W-weight in |th| > t:  frac(t) <= pi sqrt(a) e^{0.51} (1-rho)^{-1} ((1-rho)/q(t))^{a+1},  q(t)^2 = (1-rho)^2 + 4 rho sin^2(t/2).
  State relations used (exact algebra, checked on all sampled states with L >= 4.45): 0.49 d <= 1 - rho = a/N <= 0.51 d,  T d^2 >= 3.9 a^2.
  P(B < 1 - b) <= E(1-B)^998 / b^998 = (a)_998/(N+s)_998 / b^998 <= (2a/N / b)^998 <= (1.02 d / b)^998.
Dressing (k = 2 only, TAIL_CERT §4 last line, used with m = 0): |tau D_b R(zeta)| <= 40 sup_{|zeta'-zeta| = (1-rho)/2} |R(zeta')|; the circles stay in
  |zeta| < 1 and move |1 - zeta| by at most delta = 0.255 d.  Hence |V_b R| <= 41 sup over the delta-enlarged regions.
Regions (xi = B rho e^{i th}), with lower/upper bounds of |1 - xi| after enlargement:
  G1: B >= 1-b1, |th| <= th0:   |1-xi| <= om + th0 + b1 + delta              (tail series)
  G2: B >= 1-b1, th0 < |th| <= th1:  <= om + th1 + b1 + delta, weight frac(th0) (tail series)
  tier k: 1-b_{k+1} <= B < 1-b_k, |th| <= th1:  b_k - delta <= |1-xi| <= om + th1 + b_{k+1} + delta,   weight (1.02d/b_k)^998   (|R| <= |F| + |T_J|)
  last tier: B < 1 - b_K:  b_K - delta <= |1-xi| <= 2 + delta,   weight (1.02 d/b_K)^998
  Th: |th| > th1:  0.49d - delta <= |1-xi| <= 2 + delta,  weight frac(th1)
  |F_{p,nu}(xi)| <= (pi/2) e^{pi|nu|/2} |1-xi|^{-p} (FULLPOS3 Lemma B);  |(1-xi)^{-+ i w}| <= e^{pi |Re w|/2} |1-xi|^{+-Im w}.
  |V_b(R)| <= 1.6 * 41 * [sup_G1 + frac(th0) sup_G2 + sum_k P_k sup_k + frac(th1) sup_Th].
Real nu >= 0.4: |B(+-nu)| <= B_abs_real(nu_lo) (rt_cert), |c_j(+-nu)| <= prod_{i<j} ((i+1/2)^2 + nu_hi^2)/(sqrt((i+1)^2 + 4 nu_lo^2)(i+1)).
Small nu in [0, 0.4] (pole of B at 0):  B(nu) = E(nu)/(2 i nu), E(nu) = Gamma(1+2i nu)/Gamma(1/2+i nu)^2 entire near 0, so with h = E Phi,
  sum_{+-} B(+-nu) Phi(+-nu) = (h(nu) - h(-nu))/(2 i nu) = (1/2i nu) int_{-nu}^{nu} h'(t) dt,  |.| <= sup_{|t|<=nu} |h'(t)| <= (1/r') sup_{Wbox} |E Phi|,
  Wbox = {|Re w| <= 0.4 + r', |Im w| <= r'}, r' = 0.05; the same c_j bound with complex w: |a0+i|^2 <= (i+1/2+|y|)^2 + x^2, |b0+i| >= i+1-2|y|.
  (Applies to the tail series on G and to T_J on the tiers / Th; F itself has no pole.)
Tail ratio: sup_{j>=J} |c_{j+1}/c_j| <= (1 + (x+|y|)/(J+1))^2 / (1 - 2|y|/(J+1)).
"""
import math, os, sys, json
from flint import arb, acb, ctx
ctx.prec = 96
HERE = os.path.dirname(os.path.abspath(__file__))
ADS4 = os.path.abspath(os.path.join(HERE, '..', '..'))
I = acb(0, 1)
FSUP = 21.33815             # sup |f| on [0, 1300] (TAIL_CERT, as used in rt_rem.py / tc_remainder.py)
RPRIME = 0.05


def logsum(xs):
    xs = [x for x in xs if x > -1e300]
    m = max(xs)
    return m + math.log(sum(math.exp(x - m) for x in xs))


def lcs(x_lo, x_hi, y, J):
    """log upper bounds of |c_j(+-w)|, j = 0..J, Re w in [x_lo, x_hi] (x_lo >= 0), |Im w| <= y (float list of logs, Arb-rounded)."""
    out = [0.0]
    c = arb(0); h2 = arb(x_hi)**2; l2 = arb(x_lo)**2; yy = arb(y)
    for i in range(J):
        num = (arb(i) + arb(0.5) + yy)**2 + h2
        if y == 0:
            den = ((arb(i + 1))**2 + 4*l2).sqrt()*(i + 1)
        else:
            den = (arb(i + 1) - 2*yy)*(i + 1)
        c = c + (num/den).log()
        out.append(float(c.upper()))
    return out


def lsum_pow(lc, J, e0, lo, hi):
    """log of sum_{j<J} |c_j| max(hi^{j-e0}, lo^{j-e0})"""
    llo, lhi = math.log(lo), math.log(hi)
    return logsum([lc[j] + max((j - e0)*llo, (j - e0)*lhi) for j in range(J)])


def lfrac(th, d, a=998, d_lo=None):
    d_lo = d if d_lo is None else d_lo
    om = 0.51*d
    q = math.sqrt(om*om + 4*(1 - om)*math.sin(th/2)**2)
    return math.log(math.pi*math.sqrt(a)) + 0.51 - math.log(0.49*d_lo) + (a + 1)*math.log(om/q)


def log_loc3(p, J, L, lc, lpref, lexp, y, b1, th0, th1, dress=True, L_hi=None):
    """log bound of |V_b(R_J)| for the frequency cell encoded by (lc, lpref, lexp, y):
    lc = log|c_j| list, lpref = log of the prefactor of the j-series part (log 2|B|, or log(sup|E|/r') in small mode),
    lexp = pi x_hi/2, y = |Im w| bound (0 for real cells).  Returns (log bound, parts)."""
    # uniform over L in [L, L_hi]: d = d_max (= e^{-L}) in every upper radius / weight (all increasing in d), d_lo = e^{-L_hi} in lower radii
    d = math.exp(-L); d_lo = math.exp(-(L if L_hi is None else L_hi)); om = 0.51*d; e0 = p - 0.5
    dl = 0.255*d if dress else 0.0
    lo_all = 0.49*d_lo - dl
    lyf = y*abs(math.log(lo_all))             # |1-xi|^{-+y} <= lo_all^{-y}  (all |1-xi| in [lo_all, 2.01])
    x_eff = (2*lexp/math.pi)
    Rt = (1 + (x_eff + y)/(J + 1))**2/(1 - 2*y/(J + 1))
    lF = lambda lo: math.log(math.pi/2) + lexp - p*math.log(lo)
    r0 = om + th0 + b1 + dl; r1 = om + th1 + b1 + dl
    if r1*Rt >= 0.97:
        return None
    base = lpref + lexp + lyf
    lG1 = base + lc[J] + (J - e0)*math.log(r0) - math.log(1 - r0*Rt)
    lG2 = base + lc[J] + (J - e0)*math.log(r1) - math.log(1 - r1*Rt) + lfrac(th0, d, d_lo=d_lo)
    # tiers
    bs = [b1]
    while bs[-1] < 0.95:
        bs.append(min(0.95, bs[-1]*1.6))
    tiers = []
    for k in range(len(bs)):
        lo = bs[k] - dl
        hi = (om + th1 + bs[k + 1] + dl) if k + 1 < len(bs) else 2 + dl
        lw = 998*math.log(1.02*d/bs[k])
        lT = base + lsum_pow(lc, J, e0, lo, hi)
        tiers.append(lw + logsum([lF(lo), lT]))
    lT = base + lsum_pow(lc, J, e0, lo_all, 2 + dl)
    lTh = logsum([lF(lo_all), lT]) + lfrac(th1, d, d_lo=d_lo)
    parts = dict(G1=lG1, G2=lG2, tiers=logsum(tiers), Th=lTh)
    tot = logsum(list(parts.values()))
    return math.log(1.6) + (math.log(41) if dress else 0.0) + tot + 1e-9, parts


PARAMS = [(b1, th0, th1) for b1 in (0.04, 0.06, 0.08, 0.11, 0.15) for th0 in (0.02, 0.04) for th1 in (0.1, 0.15, 0.2, 0.3)]


def best_bound(p, J, L, lc, lpref, lexp, y, dress, L_hi=None):
    best = None
    for (b1, th0, th1) in PARAMS:
        r = log_loc3(p, J, L, lc, lpref, lexp, y, b1, th0, th1, dress, L_hi)
        if r is not None and (best is None or r[0] < best[0]):
            best = (r[0], (b1, th0, th1), r[1])
    return best


_Esup = None
def Esup_small():
    """sup |E(w)|, E = Gamma(1+2iw)/Gamma(1/2+iw)^2, over |Re w| <= 0.45, |Im w| <= 0.05 (Arb boxes)."""
    global _Esup
    if _Esup is None:
        m = 0.0
        nx = 36
        for k in range(nx):
            xm = -0.45 + (k + 0.5)*0.9/nx
            for ym in (-0.025, 0.025):
                w = acb(arb(xm, 0.45/nx), arb(ym, 0.025))
                E = (1 + 2*I*w).gamma()/((acb(0.5) + I*w).gamma())**2
                m = max(m, float(abs(E).upper()))
        _Esup = m
    return _Esup


def cell_bound(p, J, L, nu_lo, nu_hi, dress, L_hi=None):
    """log bound of |V_b(R_J)(nu)| uniformly for nu in [nu_lo, nu_hi]; small mode if nu_hi <= 0.4."""
    if nu_hi <= 0.4 + 1e-12:
        y = RPRIME
        lc = lcs(0.0, 0.4 + RPRIME, y, J + 1)
        lpref = math.log(Esup_small()/RPRIME)
        lexp = math.pi*(0.4 + RPRIME)/2
    else:
        import rt_cert as RC
        y = 0.0
        lc = lcs(nu_lo, nu_hi, 0.0, J + 1)
        lpref = math.log(2*RC.B_abs_real(nu_lo))
        lexp = math.pi*nu_hi/2
    return best_bound(p, J, L, lc, lpref, lexp, y, dress, L_hi)


def tail_mass(T):
    sys.path.insert(0, os.path.join(ADS4, 'code', 'contour'))
    import cb_common as KC, cb_lines as CL
    global _W
    try:
        _W
    except NameError:
        _W = KC.load()
    v, _ = CL.tail(_W, 0.0, T)
    return float(v.log().upper())


def integrals(L, J, verbose=False, L_hi=None):
    """lI2 = log int_0^oo |f| |V_b R_J| (k = 2, p = 3, dressed);  lI3 = log int_0^1 |L(R'_J)|(Lam p) dp (k = 3, p = 4, undressed).
    |f| <= FSUP on [0, 1030]; for nu >= 1030 the mass of [T, T+1] is <= tail(T) (CONTOUR); beyond 1300 the Gaussian tail
    (tail(1300) = e^{-45450}) against the crude growth of the bound (checked: bound(2000) - bound(1300) < 8000)."""
    LAM = math.sqrt(998*999)
    ls2 = []; worst2 = (-1e9, None)
    b = cell_bound(3, J, L, 0.0, 0.4, True, L_hi)
    ls2.append(math.log(0.4*FSUP) + b[0]); worst2 = max(worst2, (b[0], 0.0))
    edges = [0.4, 1.0] + [float(k) for k in range(2, 1031)]
    for n1, n2 in zip(edges[:-1], edges[1:]):
        b = cell_bound(3, J, L, n1, n2, True, L_hi)
        ls2.append(math.log(FSUP*(n2 - n1)) + b[0])
        if b[0] > worst2[0]: worst2 = (b[0], n1, b[1])
    for T in range(1030, 1300):
        b = cell_bound(3, J, L, float(T), float(T + 1), True, L_hi)
        ls2.append(tail_mass(T) + b[0])
    b1300 = cell_bound(3, J, L, 1300.0, 1301.0, True, L_hi)[0]
    b2000 = cell_bound(3, J, L, 2000.0, 2001.0, True, L_hi)[0]
    assert b2000 - b1300 < 8000 and tail_mass(1300) + b2000 < -30000
    lI2 = logsum(ls2)
    # k = 3: nu = Lam p, p in [0, 1]
    ls3 = []; worst3 = (-1e9, None)
    b = cell_bound(4, J, L, 0.0, 0.4, False, L_hi)
    ls3.append(math.log(0.4/LAM) + b[0])
    k = 0
    lo = 0.4
    while lo < LAM:
        hi = min(LAM, math.floor(lo) + 1.0)
        b = cell_bound(4, J, L, lo, hi, False, L_hi)
        ls3.append(math.log((hi - lo)/LAM) + b[0])
        if b[0] > worst3[0]: worst3 = (b[0], lo, b[1])
        lo = hi
    lI3 = logsum(ls3)
    if verbose:
        print(f'  L={L} J={J}: lI2 = {lI2:.1f} (worst cell {worst2}), lI3 = {lI3:.1f} (worst {worst3})', flush=True)
    return lI2, lI3


if __name__ == '__main__':
    sys.path.insert(0, HERE)
    L = float(sys.argv[1]); J = int(sys.argv[2])
    for nu in (0.0, 1.0, 500.0, 1000.0, 1029.0, 1040.0, 1146.0):
        if nu == 0.0:
            b = cell_bound(3, J, L, 0.0, 0.4, True)
        else:
            b = cell_bound(3, J, L, nu, nu + 1, True)
        print(f'L={L} J={J} nu={nu}: log|V_b R| <= {b[0]:.1f} params {b[1]} parts ' + ', '.join(f'{k}: {v:.0f}' for k, v in b[2].items()), flush=True)
