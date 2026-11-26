"""
tailcert/tc_remainder.py (TAIL_CERT.md §4): the connection-series remainder R_{J0}, by the localisation lemma (Lemma Loc).
  F_{p,nu} = (1-z)^{-p} P_{-1/2+i nu}((1+z)/(1-z)) = sum_{+-} B(+-nu) sum_{j<J0} c_j(+-nu)(1-z)^{j-e0 -+ i nu} + R_{J0},  e0 = p - 1/2.
Lemma Loc [P, TAIL_CERT §4.1]: if |R(xi)| <= K |1-xi|^m on |xi| <= rho2 (rho = n/(n+a-m), rho2 = (1+rho)/2), 2 <= m <= a - 3, then
  |V_b(R)| <= 2.53 K (1.0006 w)^{m-1} Q_m (1 + 40*1.5^m),  Q_m = ((998+m+1)/(998-m-1))^{m+1}     (all tail states, w <= e^{-L0}).
K(nu) [P]: real |nu| >= 0.4:
  near (|1-xi| <= 1/2): |R| <= 2|B| e^{pi nu/2} |1-xi|^{J0-e0} sum_i |c_{J0+i}| 2^{-i}   (branch series tail; |arg(1-xi)| < pi/2);
  far  (|1-xi| > 1/2): |F_{p,nu}| <= 2^p (pi/2) e^{pi nu/2} (FULLPOS3 Lemma B), |T_{J0}| <= 2|B| e^{pi nu/2} sum_{j<J0}|c_j| 2^{|j-e0|};
  K = max(near * 2^{-(J0-e0-m)}, far * 2^m).
  |nu| <= 0.4: maximum principle in nu on |nu| = 0.4 (R analytic there: the nu = 0 poles of B(+-nu) cancel, c_j(+-nu) poles at |nu| >= 1/2);
  |F_{p,nu}(xi)| <= |1-xi|^{-p} 2^{0.9} e^{0.2 pi}/(1-|xi|) (Pfaff form (1-xi)^{1/2-i nu} 2F1(1/2-i nu,1/2-i nu;1;xi), |(0.9)_k/k!| <= 1),
  1/(1-rho2) <= 4.0002/((1-m/998) w);  branch terms with |(1-xi)^{-+i nu}| <= |1-xi|^{-0.4} e^{0.2 pi} (|1-xi| <= 1), 2^{0.4} e^{0.2 pi} (else).
Weights: k = 2: |f(nu)| <= sup|psi| + sup|r| (EPSW) for nu <= 1300, CONTOUR tail formula (Y = 0) beyond; contribution relative to the
  main term: w^{5/2} int |f| |V(R)| dnu / I_W(L).   k = 3: Lam^2 int |psi3(p)| |L(R)| dp, prefactor w^3 tau, see tc_k3.py.
Usage: python3 tc_remainder.py L0 J0
"""
import os, sys, json, math
from fractions import Fraction as Fr
from flint import arb, acb, ctx
ctx.prec = 128
HERE = os.path.dirname(os.path.abspath(__file__))
ADS4 = os.path.abspath(os.path.join(HERE, '..', '..'))
PI = arb.pi()


def up(x):
    return float(arb(x).upper())


def lup(x):
    """log upper bound of a positive arb/float."""
    return float(arb(arb(x).upper()).log().upper())


def B_real_sup(nu1):
    """sup_{nu >= nu1} |B(nu)| for real nu1 > 0: |B(nu)|^2 = coth(pi nu)/(4 pi nu) (|Gamma(2i nu)|^2 = pi/(2 nu sinh 2 pi nu),
    |Gamma(1/2+i nu)|^2 = pi/cosh(pi nu)), decreasing in nu."""
    n = arb(nu1)
    return ((PI*n).coth()/(4*PI*n)).sqrt()


def c_seq(nu, J):
    """|c_j(nu)| upper bounds, j = 0..J (nu acb ball)."""
    out = [arb(1)]
    c = acb(1); a = acb(0.5) - acb(0, 1)*nu; b = acb(1) - 2*acb(0, 1)*nu
    for j in range(J):
        c = c*(a + j)*(a + j)/((b + j)*(j + 1))
        out.append(arb(abs(c).upper()))
    return out


def near_tail(nu, J0, nu_abs_max):
    """sum_{i>=0} |c_{J0+i}(nu)| 2^{-i}: explicit terms until the ratio bound r_j = |a+j|^2/(|b+j|(j+1)) * 1/2 <= 0.75, then geometric."""
    c = acb(1); a = acb(0.5) - acb(0, 1)*nu; b = acb(1) - 2*acb(0, 1)*nu
    for j in range(J0):
        c = c*(a + j)*(a + j)/((b + j)*(j + 1))
    tot = arb(0); cj = arb(abs(c).upper()); j = J0; i = 0
    while True:
        tot += cj*arb(2)**(-i)
        r = abs((a + j)*(a + j)/((b + j)*(j + 1)))/2
        if r.upper() <= 0.75 and j > 2*nu_abs_max + 10:
            # remaining sum_{i'>i} <= cj 2^{-i} * sum (0.75)^k ... (ratio of consecutive terms of the weighted series <= 0.75 from here on:
            # |a+j|^2/(|b+j|(j+1)) is decreasing in j for j > 2|nu|+10)
            tot += cj*arb(2)**(-i)*3
            break
        cj = arb((cj*abs((a + j)*(a + j)/((b + j)*(j + 1)))).upper()); j += 1; i += 1
        assert i < 20000
    return tot


def B_abs(nu):
    return abs((2*acb(0, 1)*nu).gamma()*((acb(0, 1)*nu + arb(1)/2).rgamma())**2)


def logK_real(nu1, nu2, p, J0, m):
    """log of K valid for real nu in [nu1, nu2], nu1 >= 0.4."""
    e0 = arb(p) - arb(1)/2
    nu = acb(arb((nu1 + nu2)/2, (nu2 - nu1)/2))
    numax = nu2
    lB = lup(B_real_sup(nu1))
    lexp = math.pi*nu2/2
    cs = c_seq(nu, J0)
    nt = near_tail(nu, J0, numax)
    lnear = math.log(2) + lB + lexp + lup(nt) - (J0 - float(e0.mid()) - m)*math.log(2)
    # far
    sT = arb(0)
    for j in range(J0):
        sT += cs[j]*arb(2)**abs(j - float(e0.mid()) + 0)
    lfar = m*math.log(2) + math.log(2**p*math.pi/2*math.exp(0) + 0) + lexp
    lfar = max(lfar, m*math.log(2) + math.log(2) + lB + lexp + lup(sT)) + math.log(2)
    return max(lnear, lfar)


def logK_small(p, J0, m, L0):
    """log K valid for |nu| <= 0.4 (max principle on |nu| = 0.4), for every w <= e^{-L0} the far F-bound uses w = e^{-L}: returns
    (lK0, coefficient of L) so that log K <= lK0 + cL * L."""
    e0 = p - 0.5
    lbest_near = -1e300; lbest_T = -1e300
    nb = 64
    for k in range(nb):
        th = 2*PI*(arb(k) + arb(0.5, 0.5))/nb
        nu = acb(0.4*th.cos(), 0.4*th.sin())
        for sgn in (1, -1):
            v = sgn*nu
            lB = lup(B_abs(v))
            cs = c_seq(v, J0)
            nt = near_tail(v, J0, 1.0)
            lbest_near = max(lbest_near, lB + lup(nt))
            sT = arb(0)
            for j in range(J0):
                sT += cs[j]*arb(2)**abs(j - e0)
            lbest_T = max(lbest_T, lB + lup(sT))
    ex = 0.2*math.pi + 0.4*math.log(2)
    lnear = math.log(2) + lbest_near + ex - (J0 - e0 - 0.4 - m)*math.log(2)
    # far: |F| <= 2^p 2^{0.9} e^{0.2pi} * 4.0002/((1-m/998) w)  -> log = const + L
    lF = p*math.log(2) + 0.9*math.log(2) + 0.2*math.pi + math.log(4.0002/(1 - m/998))
    lfar_T = m*math.log(2) + math.log(2) + lbest_T + ex
    lK0 = max(lnear, m*math.log(2) + lF + math.log(2), lfar_T) + math.log(2)
    return lK0, 1.0          # log K <= lK0 + 1.0 * L  (only the F-part carries the L; taking it for all parts is an upper bound)


def log_loc_factor(m, L):
    """log of 2.53 (1.0006 w)^{m-1} Q_m (1 + 40*1.5^m)."""
    Q = (m + 1)*math.log((998 + m + 1)/(998 - m - 1))
    return math.log(2.53) + (m - 1)*(math.log(1.0006) - L) + Q + math.log(1 + 40*1.5**m)


def run_k2(L0, J0):
    m = J0 - 4
    out = []
    def pr(s):
        print(s, flush=True); out.append(s)
    # |f| bound on [0, 1300]
    nf = json.load(open(os.path.join(ADS4, 'prl', 'data', 'nf02_psi.json')))
    K0c = Fr('27.48509817')
    p2 = [Fr(c)/K0c for c in nf['psi2_exact']]
    psimax = arb(0)
    for k in range(400):
        pb = arb((k + 0.5)/400, 0.5/400)
        v = abs(sum((arb(c.numerator)/c.denominator*pb**i for i, c in enumerate(p2)), arb(0)))
        psimax = v if v > psimax else psimax
    ew = json.load(open(os.path.join(ADS4, 'code', 'epswcert', 'cache', 'ew_cert_1000_re_rep.json')))
    Fsup = up(psimax) + ew['supr']*1.0001
    pr(f'k=2: sup|f| on the real axis <= sup|psi| + sup|r| = {Fsup:.5f};  J0 = {J0}, m = {m}')
    # small nu
    lK0, cL = logK_small(3, J0, m, L0)
    # int_0^0.4 |f| K <= 0.4 Fsup K
    terms = [(math.log(0.4*Fsup) + lK0, cL)]
    # [0.4, 1300] in cells of width 1 (0.4..1, then unit cells)
    edges = [0.4, 1.0] + [float(k) for k in range(2, 1301)]
    lmax = -1e300; lsum_list = []
    for n1, n2 in zip(edges[:-1], edges[1:]):
        lK = logK_real(n1, n2, 3, J0, m)
        lsum_list.append(math.log(Fsup*(n2 - n1)) + lK)
    # beyond 1300: CONTOUR tail formula (Y = 0) per unit cell up to 6000; e^{pi nu/2 + ...} K growth vs Gaussian decay
    sys.path.insert(0, os.path.join(ADS4, 'code', 'contour'))
    import cb_common as KC, cb_lines as CL
    W = KC.load()
    for T in range(1300, 6000, 10):
        intf, _ = CL.tail(W, 0.0, T)
        lK = logK_real(T, T + 10, 3, J0, m)
        lsum_list.append(lup(intf) + lK)
    mx = max(lsum_list)
    lI = mx + math.log(sum(math.exp(x - mx) for x in lsum_list))
    pr(f'k=2: log int_[0.4,6000] |f| K <= {lI:.1f};  small-nu part log <= {lK0 + math.log(0.4*Fsup):.1f} + {cL} L')
    # beyond 6000: Gaussian factor e^{-(T-1029)^2/1.62} vs K <= e^{(pi/2 + 2) T}: negligible (checked at T = 6000 in the log)
    W_I = json.load(open(os.path.join(ADS4, 'code', 'contour', 'cache', 'cb_w.json')))['cells']
    rows = []
    for L in [L0 + 0.5*k for k in range(0, 200)]:
        IW = sum(mn*(math.exp(-L*y0) - math.exp(-L*y1))/L for y0, y1, mn, mx_ in W_I)
        lrel = 2.5*(-L) + log_loc_factor(m, L) + max(lI, lK0 + cL*L + math.log(0.4*Fsup)) + math.log(2) - math.log(IW)
        rows.append((L, lrel))
    pr(f'k=2 remainder / main: log10 <= {max(r[1] for r in rows)/math.log(10):.1f} on L in [{L0}, {L0+99.5}] (decreasing in L)')
    return out, rows


if __name__ == '__main__':
    L0 = float(sys.argv[1]); J0 = int(sys.argv[2])
    out, rows = run_k2(L0, J0)
    open(os.path.join(HERE, 'cache', f'tc_remainder_k2_L{L0}_J{J0}.txt'), 'w').write('\n'.join(out) + '\n')
